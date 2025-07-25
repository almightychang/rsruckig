//! Calculation of a state-to-state trajectory.
use crate::error::{RuckigError, RuckigErrorHandler};
use crate::util::DataArrayOrVec;
use crate::{
    block::Block,
    input_parameter::{ControlInterface, DurationDiscretization, InputParameter, Synchronization},
    position_first_step1::PositionFirstOrderStep1,
    position_first_step2::PositionFirstOrderStep2,
    position_second_step1::PositionSecondOrderStep1,
    position_second_step2::PositionSecondOrderStep2,
    position_third_step1::PositionThirdOrderStep1,
    position_third_step2::PositionThirdOrderStep2,
    profile::{ControlSigns, Direction, Profile, ReachedLimits},
    result::RuckigResult,
    trajectory::Trajectory,
    velocity_second_step1::VelocitySecondOrderStep1,
    velocity_second_step2::VelocitySecondOrderStep2,
    velocity_third_step1::VelocityThirdOrderStep1,
    velocity_third_step2::VelocityThirdOrderStep2,
};

#[derive(Debug)]
pub struct TargetCalculator<const DOF: usize> {
    eps: f64,
    return_error_at_maximal_duration: bool,
    new_phase_control: DataArrayOrVec<f64, DOF>,
    pd: DataArrayOrVec<f64, DOF>,
    possible_t_syncs: Vec<f64>,
    idx: Vec<usize>,
    blocks: DataArrayOrVec<Block, DOF>,
    inp_min_velocity: DataArrayOrVec<f64, DOF>,
    inp_min_acceleration: DataArrayOrVec<f64, DOF>,
    inp_per_dof_control_interface: DataArrayOrVec<ControlInterface, DOF>,
    inp_per_dof_synchronization: DataArrayOrVec<Synchronization, DOF>,
    pub degrees_of_freedom: usize,
}

impl<const DOF: usize> TargetCalculator<DOF> {
    pub fn new(dofs: Option<usize>) -> Self {
        Self {
            blocks: DataArrayOrVec::new(dofs, Block::default()),
            inp_min_velocity: DataArrayOrVec::new(dofs, 0.0),
            inp_min_acceleration: DataArrayOrVec::new(dofs, 0.0),
            inp_per_dof_control_interface: DataArrayOrVec::new(dofs, ControlInterface::default()),
            inp_per_dof_synchronization: DataArrayOrVec::new(dofs, Synchronization::default()),
            new_phase_control: DataArrayOrVec::new(dofs, 0.0),
            pd: DataArrayOrVec::new(dofs, 0.0),
            possible_t_syncs: vec![0.0; 3 * dofs.unwrap_or(DOF) + 1],
            idx: vec![0; 3 * dofs.unwrap_or(DOF) + 1],
            eps: f64::EPSILON,
            return_error_at_maximal_duration: true,
            degrees_of_freedom: dofs.unwrap_or(DOF),
        }
    }

    // Allowing mutable reference to self for the sake of better performance.
    #[allow(clippy::wrong_self_convention)]
    fn is_input_collinear(
        &mut self,
        inp: &InputParameter<DOF>,
        limiting_direction: Direction,
        limiting_dof: usize,
    ) -> bool {
        // Check that vectors pd, v0, a0, vf, af are collinear
        for dof in 0..self.degrees_of_freedom {
            self.pd[dof] = inp.target_position[dof] - inp.current_position[dof];
        }

        let mut scale_vector: Option<&DataArrayOrVec<f64, DOF>> = None;
        let mut scale_dof: Option<usize> = None;
        for dof in 0..self.degrees_of_freedom {
            if self.inp_per_dof_synchronization[dof] != Synchronization::Phase {
                continue;
            }

            if self.inp_per_dof_control_interface[dof] == ControlInterface::Position
                && self.pd[dof].abs() > self.eps
            {
                scale_vector = Some(&self.pd);
                scale_dof = Some(dof);
                break;
            } else if inp.current_velocity[dof].abs() > self.eps {
                scale_vector = Some(&inp.current_velocity);
                scale_dof = Some(dof);
                break;
            } else if inp.current_acceleration[dof].abs() > self.eps {
                scale_vector = Some(&inp.current_acceleration);
                scale_dof = Some(dof);
                break;
            } else if inp.target_velocity[dof].abs() > self.eps {
                scale_vector = Some(&inp.target_velocity);
                scale_dof = Some(dof);
                break;
            } else if inp.target_acceleration[dof].abs() > self.eps {
                scale_vector = Some(&inp.target_acceleration);
                scale_dof = Some(dof);
                break;
            }
        }

        if scale_dof.is_none() {
            return false; // Zero everywhere is in theory collinear, but that trivial case is better handled elsewhere
        }

        let scale = scale_vector.unwrap()[scale_dof.unwrap()];
        let pd_scale = self.pd[scale_dof.unwrap()] / scale;
        let v0_scale = inp.current_velocity[scale_dof.unwrap()] / scale;
        let vf_scale = inp.target_velocity[scale_dof.unwrap()] / scale;
        let a0_scale = inp.current_acceleration[scale_dof.unwrap()] / scale;
        let af_scale = inp.target_acceleration[scale_dof.unwrap()] / scale;

        let scale_limiting = scale_vector.unwrap()[limiting_dof];
        let mut control_limiting = if limiting_direction == Direction::UP {
            inp.max_jerk[limiting_dof]
        } else {
            -inp.max_jerk[limiting_dof]
        };
        if inp.max_jerk[limiting_dof].is_infinite() {
            control_limiting = if limiting_direction == Direction::UP {
                inp.max_acceleration[limiting_dof]
            } else {
                self.inp_min_acceleration[limiting_dof]
            };
        }

        for dof in 0..self.degrees_of_freedom {
            if self.inp_per_dof_synchronization[dof] != Synchronization::Phase {
                continue;
            }

            let current_scale = scale_vector.unwrap()[dof];
            if (self.inp_per_dof_control_interface[dof] == ControlInterface::Position
                && (self.pd[dof] - pd_scale * current_scale).abs() > self.eps)
                || (inp.current_velocity[dof] - v0_scale * current_scale).abs() > self.eps
                || (inp.current_acceleration[dof] - a0_scale * current_scale).abs() > self.eps
                || (inp.target_velocity[dof] - vf_scale * current_scale).abs() > self.eps
                || (inp.target_acceleration[dof] - af_scale * current_scale).abs() > self.eps
            {
                return false;
            }

            self.new_phase_control[dof] = control_limiting * current_scale / scale_limiting;
        }

        true
    }

    fn synchronize(
        &mut self,
        t_min: Option<f64>,
        t_sync: &mut f64,
        limiting_dof: &mut Option<usize>,
        profiles: &mut DataArrayOrVec<Profile, { DOF }>,
        discrete_duration: bool,
        delta_time: f64,
    ) -> bool {
        // Check for (degrees_of_freedom == 1 && !t_min && !discrete_duration) is now outside

        // Possible t_syncs are the start times of the intervals and optional t_min
        let mut any_interval = false;
        for dof in 0..self.degrees_of_freedom {
            // Ignore DoFs without synchronization here
            if self.inp_per_dof_synchronization[dof] == Synchronization::None {
                self.possible_t_syncs[dof] = 0.0;
                self.possible_t_syncs[self.degrees_of_freedom + dof] = f64::INFINITY;
                self.possible_t_syncs[2 * self.degrees_of_freedom + dof] = f64::INFINITY;
                continue;
            }

            self.possible_t_syncs[dof] = self.blocks[dof].t_min;
            self.possible_t_syncs[self.degrees_of_freedom + dof] =
                if let Some(a) = &self.blocks[dof].a {
                    a.right
                } else {
                    f64::INFINITY
                };
            self.possible_t_syncs[2 * self.degrees_of_freedom + dof] =
                if let Some(b) = &self.blocks[dof].b {
                    b.right
                } else {
                    f64::INFINITY
                };
            any_interval |= self.blocks[dof].a.is_some() || self.blocks[dof].b.is_some();
        }
        self.possible_t_syncs[3 * self.degrees_of_freedom] = t_min.unwrap_or(f64::INFINITY);
        any_interval |= t_min.is_some();

        if discrete_duration {
            for possible_t_sync in &mut self.possible_t_syncs {
                if possible_t_sync.is_infinite() {
                    continue;
                }

                let remainder = *possible_t_sync % delta_time; // in [0, delta_time)
                if remainder > self.eps {
                    *possible_t_sync += delta_time - remainder;
                }
            }
        }

        // Test them in sorted order
        // Setting up the range for `idx_end`
        let idx_end = if any_interval {
            self.idx.len()
        } else {
            self.degrees_of_freedom
        };

        // Initialize the range similar to `std::iota`
        for i in 0..idx_end {
            self.idx[i] = i;
        }

        // Sort the values in the range
        self.idx[0..idx_end].sort_by(|&i, &j| {
            self.possible_t_syncs[i]
                .partial_cmp(&self.possible_t_syncs[j])
                .unwrap()
        });

        // Start at last tmin (or worse)
        for &i in &self.idx[(self.degrees_of_freedom - 1)..] {
            let possible_t_sync = self.possible_t_syncs[i];
            let mut is_blocked = false;
            for dof in 0..self.degrees_of_freedom {
                if self.inp_per_dof_synchronization[dof] == Synchronization::None {
                    continue; // inner dof loop
                }
                if self.blocks[dof].is_blocked(possible_t_sync) {
                    is_blocked = true;
                    break; // inner dof loop
                }
            }
            if is_blocked || possible_t_sync < t_min.unwrap_or(0.0) || possible_t_sync.is_infinite()
            {
                continue;
            }

            *t_sync = possible_t_sync;
            if i == 3 * self.degrees_of_freedom {
                // Optional t_min
                *limiting_dof = None;
                return true;
            }

            let div = i / self.degrees_of_freedom;
            *limiting_dof = Some(i % self.degrees_of_freedom);
            match div {
                0 => {
                    profiles[limiting_dof.unwrap()] =
                        self.blocks[limiting_dof.unwrap()].p_min.clone();
                }
                1 => {
                    profiles[limiting_dof.unwrap()] = self.blocks[limiting_dof.unwrap()]
                        .a
                        .clone()
                        .unwrap()
                        .profile;
                }
                2 => {
                    profiles[limiting_dof.unwrap()] = self.blocks[limiting_dof.unwrap()]
                        .b
                        .clone()
                        .unwrap()
                        .profile;
                }
                _ => {}
            }
            return true;
        }

        false
    }

    /// Calculate the time-optimal waypoint-based trajectory.
    pub fn calculate<T: RuckigErrorHandler>(
        &mut self,
        inp: &InputParameter<DOF>,
        traj: &mut Trajectory<DOF>,
        delta_time: f64,
    ) -> Result<RuckigResult, RuckigError> {
        for dof in 0..self.degrees_of_freedom {
            let p = &mut traj.profiles[0][dof];

            self.inp_min_velocity[dof] = inp
                .min_velocity
                .as_ref()
                .map_or(-inp.max_velocity[dof], |v| v[dof]);

            self.inp_min_acceleration[dof] = inp
                .min_acceleration
                .as_ref()
                .map_or(-inp.max_acceleration[dof], |v| v[dof]);

            self.inp_per_dof_control_interface =
                DataArrayOrVec::new(Some(self.degrees_of_freedom), inp.control_interface.clone());
            if let Some(per_dof_control_interface) = &inp.per_dof_control_interface {
                for (dof, value) in per_dof_control_interface.iter().enumerate() {
                    *self.inp_per_dof_control_interface.get_mut(dof).unwrap() = value.clone();
                }
            }

            self.inp_per_dof_synchronization =
                DataArrayOrVec::new(Some(self.degrees_of_freedom), inp.synchronization.clone());
            if let Some(per_dof_synchronization) = &inp.per_dof_synchronization {
                for (dof, value) in per_dof_synchronization.iter().enumerate() {
                    *self.inp_per_dof_synchronization.get_mut(dof).unwrap() = value.clone();
                }
            }

            if !inp.enabled[dof] {
                if let Some(last) = p.p.last_mut() {
                    *last = inp.current_position[dof];
                }
                if let Some(last) = p.v.last_mut() {
                    *last = inp.current_velocity[dof];
                }
                if let Some(last) = p.a.last_mut() {
                    *last = inp.current_acceleration[dof];
                }
                if let Some(last) = p.t_sum.last_mut() {
                    *last = 0.0;
                }

                self.blocks[dof].t_min = 0.0;
                self.blocks[dof].a = None;
                self.blocks[dof].b = None;
                continue;
            }

            // Calculate brake (if input exceeds or will exceed limits)
            p.brake.get_position_brake_trajectory(
                inp.current_velocity[dof],
                inp.current_acceleration[dof],
                inp.max_velocity[dof],
                inp.min_velocity
                    .as_ref()
                    .and_then(|v| v.get(dof))
                    .cloned()
                    .unwrap_or(-inp.max_velocity[dof]),
                inp.max_acceleration[dof],
                inp.min_acceleration
                    .as_ref()
                    .and_then(|v| v.get(dof))
                    .cloned()
                    .unwrap_or(-inp.max_acceleration[dof]),
                inp.max_jerk[dof],
            );

            // Finalize pre & post-trajectories
            p.brake.finalize(&mut p.p[0], &mut p.v[0], &mut p.a[0]);

            let mut step1 = PositionThirdOrderStep1::new(
                p.p[0],
                p.v[0],
                p.a[0],
                p.pf,
                p.vf,
                p.af,
                inp.max_velocity[dof],
                inp.min_velocity
                    .as_ref()
                    .map_or(-inp.max_velocity[dof], |v| v[dof]),
                inp.max_acceleration[dof],
                inp.min_acceleration
                    .as_ref()
                    .map_or(-inp.max_acceleration[dof], |v| v[dof]),
                inp.max_jerk[dof],
            );
            let found_profile = step1.get_profile(p, &mut self.blocks[dof]);

            if !found_profile {
                let has_zero_limits = inp.max_acceleration[dof] == 0.0
                    || inp
                        .min_acceleration
                        .as_ref()
                        .map_or(-inp.max_acceleration[dof], |v| v[dof])
                        == 0.0
                    || inp.max_jerk[dof] == 0.0;
                if has_zero_limits {
                    T::handle_calculator_error(&format!(
                        "zero limits conflict in step 1, dof: {} input: {}",
                        dof, inp
                    ))?;
                    return Ok(RuckigResult::ErrorZeroLimits);
                }
                T::handle_calculator_error(&format!(
                    "error in step 1, dof: {} input: {}",
                    dof, inp
                ))?;
                return Ok(RuckigResult::ErrorExecutionTimeCalculation);
            }

            traj.independent_min_durations[dof] = self.blocks[dof].t_min;
        }

        traj.duration = self.blocks[0].t_min;
        traj.profiles[0][0] = self.blocks[0].p_min.clone();
        traj.cumulative_times[0] = traj.duration;

        Ok(RuckigResult::Working)
    }
}
