# MTR 5100 Final Project — Camera-Guided RCM Intravitreal Injection Robot

**Author:** Mohamed Abdelaziz Hamdan (481019296)
**Software:** MATLAB R2024a, base MATLAB only. No toolboxes are required.

## System overview
| Part | Joints | Job |
|---|---|---|
| Stereo camera (2 calibrated pinhole cameras on a boom) | – | Locates the eye: limbus points → triangulation → circle fit → eye centre and pars plana entry point |
| 4-DOF positioning arm | q1 base yaw, q2 shoulder, q3 elbow, q4 wrist pitch | Carries the injection unit from its start pose to the eye along a collision-free RRT-Connect path (the patient's face is the obstacle), docks it, and brings it back |
| 3-DOF RCM injection robot | q1 yaw, q2 pitch (R-R), q3 insertion (P) | Inserts the needle about the remote centre of motion and injects. J2 is built as a curved guide rail centred on the RCM (arc-guide), so no bearing sits on the pitch axis next to the face; a numerical check of the robot against the face model gives a clearance of at least 1.4 mm through the whole demo |

## How to run
1. Open MATLAB and `cd` into this `code` folder.
   - **Explore the room without running the procedure:** `explore_scene`. It opens the operating room with preset views (overview, face and eye, eye close-up, RCM unit, needle entry, arm, camera stand, surgeon's view, the stereo camera's own view), free navigation (mouse: orbit / pan / zoom; keys: W/A/S/D/Q/E walk, arrows look around, R reset), robot poses and joint sliders for both robots with a live clearance read-out, a *Plan path to the eye* button (RRT-Connect from the current arm pose), and buttons that start the live demo, record the video or run the full analysis.
2. Run `main_final_project`.
   - First it plays the **live 3-D prototype demo** (about 1.5 minutes, 15 steps):
     1. The stereo camera locates the eye.
     2. The 4-DOF arm solves the docking pose with closed-form inverse kinematics; a ghost arm previews the solution.
     3. The arm carries the unit along a path planned by RRT-Connect in its 4-D joint space (patient face SDF, body, table and camera stand as obstacles, plus self-collision of the arm and the unit); the speed is limited so the protective-stop distance always fits in the clearance (speed and separation monitoring).
     4. The arm docks along the insertion axis onto the **measured** entry point.
     5. The 3-DOF RCM robot inserts the needle, injects at two sites (IOP feedback), re-orients with A*, and survives an eye twitch.
     6. The arm undocks and returns to the exact start pose along a second planned path.
   - Then it runs the full engineering analysis. Set `RUN_DEMO = false` at the top of the file to skip the demo.
3. To run only the demo: `demo_live`. Other modes:
   - `demo_live('record')`: also saves `../results/demo_procedure.mp4`, a backup video for the presentation
   - `demo_live('snapshots')`: also saves one image per step (`../figures/D01 … D15`)
   - `demo_live('fast')`: plays without real-time pacing
4. Outputs:
   - `../figures/F01 … F14*.png`: every figure in the report (200 dpi)
   - `../results/run_log.txt`: the full console log with all tables
   - `../results/results.mat`: the summary metrics and key result arrays

The analysis takes about 11 to 26 minutes, depending on the machine. All random numbers come from fixed seeds (`P.seed`), so every run gives the same results.

Batch alternative: `matlab -batch "main_final_project"`.

## DH tables
**3-DOF RCM injection robot** (base frame {0} at the entry point, x0 toward the eye centre)

| i | θi | di | ai | αi |
|---|---|---|---|---|
| 1 | q1 | 0 | 0 | −90° |
| 2 | q2 | 0 | 0 | +90° |
| 3 | 0 | q3 | 0 | 0 |

Needle direction is u = [c1 s2, s1 s2, c2], and the tip is at p = q3 · u.

**4-DOF positioning arm** (metres)

| i | θi | di | ai | αi |
|---|---|---|---|---|
| 1 | q1 | 0.400 | 0 | +90° |
| 2 | q2 | 0 | 0.300 | 0 |
| 3 | q3 | 0 | 0.260 | 0 |
| 4 | q4 | 0 | 0.120 | 0 |

The task space is the flange position plus its pitch angle φ = q2 + q3 + q4. The IK has 4 solutions (front/back × elbow up/down). The RCM unit is fixed to the flange through `P.arm.T_40`, so T_E0 = T_E4(qa) · T_40.

## File map
| File | Purpose | Course module |
|---|---|---|
| `main_final_project.m` | Main script, sections S1–S14 | all |
| `demo_live.m` | Live 3-D prototype animation of the complete procedure (15 steps) | all |
| `explore_scene.m` | Interactive walk-through of the operating room: preset views, navigation, poses, joint sliders with live clearance, on-demand arm planning, demo / analysis launch buttons | all |
| `fp_arm_plan.m` | Arm transport planner: RRT-Connect (best of 3) → shortcuts → Gaussian corner rounding → time law with speed and separation monitoring → full-mesh verification | M7, M8, M9 |
| `fp_rrt_connect.m` | Bidirectional RRT (RRT-Connect) in an n-D configuration space | M8 |
| `fp_arm_collision_model.m`, `fp_arm_clearance.m` | Point-cloud model of the arm and RCM unit taken from the drawn 3-D geometry; smallest distance to the obstacles for any arm / RCM joint values | M8 |
| `fp_arm_selfclear.m` | Self-collision clearance of the arm: links and column as capsules, the RCM unit as points; all non-adjacent pairs (unit vs column / links, link vs link) | M8 |
| `fp_arm_env.m` | Obstacle distance field: patient head (fine + coarse grids, a guaranteed lower bound of the true distance), body, table, camera stand | M8, M9 |
| `fp_eye_register.m` | Stereo camera model, DLT triangulation, plane and circle fit, entry-point estimate | M3, M4 |
| `fp_arm_fk.m`, `fp_arm_ik.m` | 4-DOF arm: DH forward kinematics and closed-form IK (4 branches, closest-branch selection) | M4, M5 |
| `fp_scene.m`, `fp_scene_update.m` | 3-D prototype model, driven by the DH frames: operating room (floor, table on a pedestal), patient, textured eye (iris, vessels, cornea) held open by a wire speculum, stereo camera on a mobile stand, 4-DOF arm (joint caps, cable hoses, E-stop), RCM unit and syringe | M1, M4, M5 |
| `fp_head_mesh.m` | Patient: the face is an anthropometric height surface (profile, nose, lips, orbits, cheeks, chin) closed by a face volume, cranium, ears and neck; meshed with `isosurface`; plus the body under the blanket. Also returns the distance fields used for collision checks | – |
| `fp_camera_image.m` | Renders the scene from the left camera's exact pose and intrinsics: the realistic camera frame shown in the demo and in F10 | M3, M4 |
| `fp_set_camera.m` | Perspective camera for the demo's cinematic shots | – |
| `fp_params.m` | All parameters and assumptions (eye, robots, camera, actuators, sensors, IOP, safety) | – |
| `fp_fk.m`, `fp_fk_geom.m` | RCM robot DH FK (with DH-error injection) and independent geometric FK | M4, M5 |
| `fp_ik.m` | RCM robot analytical IK (spherical coordinates, 4 candidates, joint limits) | M5 |
| `fp_jacobian.m` | Geometric Jacobian plus the RCM-point velocity Jacobian | M6 |
| `fp_arm_dynamics.m` | 4-DOF arm dynamics M(q), c, G with the carried RCM unit (Lagrange / Christoffel), plus reflected rotor inertia | M2, M6 |
| `fp_dynamics.m` | M(q), c(q,qd), G(q) via the Lagrange/Christoffel method, plus reflected rotor inertia | M2, M6 |
| `fp_traj.m` | Quintic and trapezoidal (synchronised) trajectories | M7 |
| `fp_simulate.m` | Closed-loop RK4 simulation (PID / computed torque, encoder quantisation, saturation, needle force) | M3, M7 |
| `fp_needle_force.m` | Sclera puncture and needle-friction force model | M3 |
| `fp_butter.m`, `fp_sosfilt.m` | Butterworth tremor filter designed by bilinear transform | M3, M7 |
| `fp_lens_check.m` | Collision/clearance checks: lens, retina, keep-out zone | M8 |
| `fp_astar.m`, `fp_rrt.m`, `fp_shortcut.m` | A*, RRT and shortcut path smoothing in C-space | M8 |
| `fp_iop_sim.m` | Injection hydraulics (Poiseuille) and IOP dynamics (Friedenwald + Goldmann) | M2, M3, M9 |
| `fp_twitch_sim.m` | Eye-twitch detection and emergency retraction | M3, M9 |

## Script sections
S1 workspace · S2 frames + camera registration (F10) · S3 FK/IK of both robots (F11) · S4 Jacobian, statics, dynamics · S5 encoder resolution · S6 tremor filter and teleoperation · S7 motion planning inside the eye · S7b arm transport planning around the patient (RRT-Connect, speed and separation monitoring, F13) · S8 full mission control · S9 actuator feasibility · S9b arm dynamics, drive sizing from a motor catalogue, power-off brakes and protective-stop braking capability (F14) · S10 RCM calibration · S11 IOP · S12 eye-twitch safety · S13 Module 9 safety case: stopping distance, protective vs emergency stop, risk register S×L before/after, residual limitations, ethics (F12) · S14 summary and FMEA

## Verification cases built into the script
- DH FK compared with the independent geometric FK, and FK(IK(p)) checked for 5000 random targets
- Arm FK(IK(T)) checked for 3000 random poses; the demo prints the docking residual and the return-to-start error
- Camera registration: Monte Carlo over 300 noise/calibration draws, compared with the true entry point
- Analytical Jacobian compared with a finite-difference Jacobian
- Dynamics: M is symmetric positive definite, G = ∂V/∂q, the passivity identity holds, and energy is conserved in free motion
- Closed-eye IOP compared with the analytical result P0·exp(K·V)
- Butterworth magnitude response at the cut-off frequency is −3 dB
- Eye-tracker false-alarm test during 20 s of steady fixation
- Patient distance field compared with the true distance to the head mesh (it must never over-estimate)
- Arm paths verified with every mesh vertex; planner tested from random start poses
- Arm dynamics: M symmetric positive definite, G = ∂V/∂q, passivity identity; the braking deceleration assumed in the M9 stopping distance is checked against the drive model
