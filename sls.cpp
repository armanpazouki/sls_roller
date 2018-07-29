#include <stdio.h>
#include <vector>
#include <cmath>
#include <string>
#include <iostream>

#include "chrono_parallel/physics/ChSystemParallel.h"
#include "chrono/utils/ChUtilsCreators.h"
#include "chrono/utils/ChUtilsGenerators.h"

//#undef CHRONO_OPENGL
#ifdef CHRONO_OPENGL
#include "chrono_opengl/ChOpenGLWindow.h"
#endif

using std::cout;
using std::endl;

using namespace chrono;
using namespace chrono::collision;

// all dimensions are in millimeters, milligrams
real container_width = 2.5;      // width of area with particles
real container_length = 25;      // length of area that roller will go over		1194mm maximum
real container_thickness = .25;  // thickness of container walls
real container_height = 2;       // height of the outer walls
real container_friction = 0;
real floor_friction = .2;
real spacer_width = 1;
real spacer_height = 1;

real roller_overlap = 1;          // amount that roller goes over the container area
real roller_length = 2.5 - .25;   // length of the roller
real roller_radius = 76.2 / 2.0;  // radius of roller
real roller_omega = 0;
real roller_velocity = -127;	  //change from -20 to -1000 with ~6 steps in between
real roller_mass = 1;
real roller_friction = .1;
real roller_cohesion = 0;
real particle_radius = 2.5 * .058 / 2.0; //.058 / 2.0;	Note: 3* = 50,000 particles; 1.5* = 500,000 particles
real particle_std_dev = .015 / 2.0;
real particle_mass = .05;
real particle_density = 0.93;
real particle_layer_thickness = 0.928; // particle_radius * 32;
real particle_friction = .52;
real rolling_friction = .1;
real spinning_friction = .1;
real gravity = -9810;  // acceleration due to gravity
                       // step size which will not allow interpenetration more than 1/6 of smallest radius
// real timestep = Abs(((particle_radius - particle_std_dev) / 3.0) / roller_velocity);
real timestep = .00005;  // step size, original = 0.00005
real time_end = 1;       // length of simulation
real current_time = 0;
int out_fps = 6000;			// original = 6000 (proportional to timestep)
int out_steps = std::ceil((1.0 / timestep) / out_fps);

int num_steps = time_end / timestep;
int max_iteration = 15;
int tolerance = 0;

int threads = 64;

std::string data_output_path = "data_sls";
std::shared_ptr<ChBody> ROLLER;
real ang = 0;

inline void RunTimeStep(ChSystemParallelNSC* mSys, const int frame) {
    auto roller_pos = ROLLER->GetPos();
    ROLLER->SetPos(ChVector<>(0, roller_radius + particle_layer_thickness + container_thickness,
                              roller_pos.z() + roller_velocity * timestep));
    ROLLER->SetPos_dt(ChVector<>(0, 0, roller_velocity));
    roller_omega = roller_velocity / roller_radius;
    ang += roller_omega * timestep;
    if (ang >= 2 * CH_C_PI) {
        ang -= 2 * CH_C_PI;
    }

    Quaternion q1;
    q1.Q_from_AngY(ang);
    Quaternion q2;
    q1 = Q_from_AngX(-ang);

    ChQuaternion<> roller_quat;
    roller_quat.Q_from_AngAxis(CH_C_PI / 2.0, ChVector<>(0, 0, 1));
    ROLLER->SetRot(q1 % roller_quat);
    ROLLER->SetWvel_loc(ChVector<>(0, roller_omega, 0));

    cout << "step " << frame << " " << ROLLER->GetPos().z() << "\n";
}

// set arguments from input
// -----------------------------------------------------------------------------
void SetArgumentsForSlsFromInput(int argc, char* argv[]) {
	int problemTypeInt = 0; // 0: SETTLING, 1: FLOWING
	if (argc > 1) {
		const char* text = argv[1];
		roller_velocity = atof(text);
	}
	if (argc > 2) {
		const char* text = argv[2];
		particle_friction = atof(text);
	}
	if (argc > 3) {
		const char* text = argv[3];
		timestep = atof(text);			// original = 0.00005
		out_fps = 6000 * (0.0005 / timestep);		// orignal = 6000
		int out_steps = std::ceil((1.0 / timestep) / out_fps);
	}
}
// -----------------------------------------------------------------------------
// jsanta35/../> sls -100 .1
int main(int argc, char* argv[]) {
	// Get problem parameters from arguments
	SetArgumentsForSlsFromInput(argc, argv);

    real roller_start = container_length + roller_radius / 3.0;

    time_end = (roller_start) / Abs(roller_velocity);

    printf("Time to run: %f %f %f\n", roller_start, time_end, timestep);

    num_steps = time_end / timestep;

	ChSystemParallelNSC* mSystem = new ChSystemParallelNSC;
    mSystem->Set_G_acc(ChVector<>(0, gravity, 0));
    mSystem->SetTimestepperType(ChTimestepper::Type::EULER_IMPLICIT_LINEARIZED);

	// **** Set number of threads.
	int max_threads = omp_get_num_procs();
	if (threads > max_threads)
		threads = max_threads;
	mSystem->SetParallelThreadNumber(threads);
	omp_set_num_threads(threads);
	cout << "Using " << threads << " threads" << endl;	// Set number of threads.
	//

    mSystem->GetSettings()->solver.tolerance = tolerance;
    mSystem->GetSettings()->solver.solver_mode = SolverMode::SPINNING;
    mSystem->GetSettings()->solver.max_iteration_normal = 30;
    mSystem->GetSettings()->solver.max_iteration_sliding = max_iteration;
    mSystem->GetSettings()->solver.max_iteration_spinning = max_iteration;
    mSystem->GetSettings()->solver.max_iteration_bilateral = 0;  // make 1000, should be about 220
    mSystem->GetSettings()->solver.compute_N = false;
    mSystem->GetSettings()->solver.alpha = 0;
    mSystem->GetSettings()->solver.cache_step_length = true;
    mSystem->GetSettings()->solver.use_full_inertia_tensor = false;
    mSystem->GetSettings()->solver.contact_recovery_speed = 180;
    mSystem->GetSettings()->solver.bilateral_clamp_speed = 1e8;
    mSystem->GetSettings()->collision.aabb_max = real3(4.50145, 77.3794, 75.8014);
    mSystem->GetSettings()->collision.aabb_min = real3(-4.50145, -0.125, -25.0014);
    mSystem->GetSettings()->collision.use_aabb_active = true;
    mSystem->ChangeSolverType(SolverType::BB);
	mSystem->SetLoggingLevel(LoggingLevel::LOG_TRACE, true);
	mSystem->SetLoggingLevel(LoggingLevel::LOG_INFO, true);

    mSystem->GetSettings()->collision.collision_envelope = particle_radius * .05;
    mSystem->GetSettings()->collision.bins_per_axis = vec3(40, 300, 400);
    mSystem->GetSettings()->collision.narrowphase_algorithm = NarrowPhaseType::NARROWPHASE_HYBRID_MPR;
    mSystem->GetSettings()->collision.fixed_bins = true;

    auto material_plate = std::make_shared<ChMaterialSurfaceNSC>();
    material_plate->SetFriction(0);
	auto PLATE = std::make_shared<ChBody>(std::make_shared<ChCollisionModelParallel>());
    utils::InitializeObject(PLATE, 100000, material_plate, ChVector<>(0, 0, 0), QUNIT, true, true, 2, 6);

    utils::AddBoxGeometry(PLATE.get(), ChVector<>(container_thickness, container_height, container_length),
                          ChVector<>(-container_width + container_thickness, container_height, 0));
    utils::AddBoxGeometry(PLATE.get(), ChVector<>(container_thickness, container_height, container_length),
                          ChVector<>(container_width - container_thickness, container_height, 0));
    utils::AddBoxGeometry(PLATE.get(), ChVector<>(container_width, container_height, container_thickness),
                          ChVector<>(0, container_height, -container_length + container_thickness));
    utils::AddBoxGeometry(PLATE.get(), ChVector<>(container_width, container_height, container_thickness),
                          ChVector<>(0, container_height, container_length - container_thickness));
    utils::AddBoxGeometry(PLATE.get(), ChVector<>(container_width, container_thickness, container_length),
                          ChVector<>(0, container_height * 2, 0));

    utils::FinalizeObject(PLATE, (ChSystemParallel*)mSystem);

    auto material_bottom = std::make_shared<ChMaterialSurfaceNSC>();
    material_bottom->SetFriction(floor_friction);
	auto BOTTOM = std::make_shared<ChBody>(std::make_shared<ChCollisionModelParallel>());
    utils::InitializeObject(BOTTOM, 100000, material_bottom, ChVector<>(0, 0, 0), QUNIT, true, true, 2, 6);
    utils::AddBoxGeometry(BOTTOM.get(), ChVector<>(container_width, container_thickness, container_length));
    utils::FinalizeObject(BOTTOM, (ChSystemParallel*)mSystem);

	ROLLER = std::make_shared<ChBody>(std::make_shared<ChCollisionModelParallel>());
    ChQuaternion<> roller_quat;
    roller_quat.Q_from_AngAxis(CH_C_PI / 2.0, ChVector<>(0, 0, 1));

    auto material_roller = std::make_shared<ChMaterialSurfaceNSC>();
    material_roller->SetFriction(roller_friction);

    utils::InitializeObject(ROLLER, 100000, material_roller,
                            ChVector<>(0, roller_radius + particle_layer_thickness + container_thickness, roller_start),
                            roller_quat, true, false, 6, 6);

    utils::AddCylinderGeometry(ROLLER.get(), roller_radius, roller_length * 2);
    utils::FinalizeObject(ROLLER, (ChSystemParallel*)mSystem);

    auto material_granular = std::make_shared<ChMaterialSurfaceNSC>();
    material_granular->SetFriction(particle_friction);
    material_granular->SetRollingFriction(rolling_friction);
    material_granular->SetSpinningFriction(spinning_friction);

    utils::Generator* gen = new utils::Generator(mSystem);

    auto m1 = gen->AddMixtureIngredient(utils::SPHERE, 1);
    m1->setDefaultSize(particle_radius);
    m1->setDefaultDensity(particle_density);
    m1->setDistributionSize(particle_radius, particle_std_dev, particle_radius - particle_std_dev,
                            particle_radius + particle_std_dev);
    m1->setDefaultMaterial(material_granular);

    gen->createObjectsBox(utils::HCP_PACK, (particle_radius + particle_std_dev) * 2, ChVector<>(0, 1.0 + particle_layer_thickness*.5, 0),
                          ChVector<>(container_width - container_thickness * 2.5, particle_layer_thickness,
                                     container_length - container_thickness * 2.5));

#ifdef CHRONO_OPENGL
    opengl::ChOpenGLWindow& gl_window = opengl::ChOpenGLWindow::getInstance();
	gl_window.Initialize(1280, 720, "Bucky", mSystem);
    gl_window.SetCamera(ChVector<>(0, 0, -10), ChVector<>(0, 0, 0), ChVector<>(0, 1, 0), 0.1);
    gl_window.Pause();
    int frame = 0;

    while (frame < num_steps) {
        if (gl_window.Active()) {
            if (gl_window.DoStepDynamics(timestep)) {
                // TimingOutput(mSystem);
                RunTimeStep(mSystem, frame);
                frame++;
            }
            gl_window.Render();
        } else {
            exit(0);
        }
    }
#else

    double time = 0, exec_time = 0;
    int sim_frame = 0, out_frame = 0, next_out_frame = 0;

    while (time < time_end) {
        mSystem->DoStepDynamics(timestep);
        if (sim_frame == next_out_frame) {
            std::cout << "write: " << out_frame << std::endl;
			utils::WriteShapesPovray(mSystem, data_output_path + "data_" + std::to_string(out_frame) + ".txt",
                                             true);
            out_frame++;
            next_out_frame += out_steps;
        }
		RunTimeStep(mSystem, sim_frame);

        // Update counters.
        time += timestep;
        sim_frame++;
        exec_time += mSystem->GetTimerStep();
    }
    cout << "==================================" << endl;
    cout << "Simulation time:   " << exec_time << endl;

#endif
}
