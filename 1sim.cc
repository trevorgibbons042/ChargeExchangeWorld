#include <iostream>

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4VisManager.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

#include "2construction.hh"
#include "3physics.hh"
#include "4action.hh"
#include "6detector.hh"
#include "10tracking.hh"
#include "G4HadronicParameters.hh"
#include "G4RunManagerFactory.hh"
#include "Randomize.hh"
#include <random>

int main(int argc, char** argv)
{ 
    auto* runManager = G4RunManagerFactory::CreateRunManager();
    //thread count comes from /run/numberOfThreads in the macro (set before /run/initialize)
    std::random_device rd;
    G4long seed = static_cast<G4long>(rd());
    G4Random::setTheSeed(seed);

    runManager->SetUserInitialization(new MyDetectorConstruction()); //construction.cc stuff
    runManager->SetUserInitialization(new MyPhysicsList()); //physics.cc stuff
    runManager->SetUserInitialization(new MyActionInitialization()); //action.cc stuff

    //initialization happens via /run/initialize in the macro / interactive branch, so /run/numberOfThreads takes effect

    G4UIExecutive *ui = nullptr;
    if(argc == 1){
        ui = new G4UIExecutive(argc, argv);
    }

    G4UImanager *UImanager = G4UImanager::GetUIpointer(); //setup for commands
    
    if(ui){
        G4VisManager *visManager = new G4VisExecutive(); //setup for visuals
        visManager -> Initialize(); //activates g4 visuals
        
        UImanager->ApplyCommand("/run/initialize");
        UImanager->ApplyCommand("/vis/open OGL"); //uses opengl
        //UImanager->ApplyCommand("/control/execute 1vis.mac");
        UImanager->ApplyCommand("/vis/drawVolume"); //draws world
        //the world is centred on the gun (z = 0) and about 68 m long, so frame only about z = 0.5 to 38.5 m (the default 37 m beamline; the scale and target values below assume zdcFrontZ = 37 m), side view with z running left to right
        UImanager->ApplyCommand("/vis/geometry/set/visibility logicWorld 0 true"); //draws the world box (off by default)
        UImanager->ApplyCommand("/vis/viewer/set/viewpointVector -1 0 0");
        UImanager->ApplyCommand("/vis/viewer/set/upVector 0 1 0");
        UImanager->ApplyCommand("/vis/viewer/set/targetPoint 0 0 19.5 m"); //centre of about 0.5 to 38.5 m, so the vertical scale just past the box end stays in frame
        UImanager->ApplyCommand("/vis/viewer/addCutawayPlane 0 0 -1 m 0 0 1"); //hides the world box and anything else behind z = -1 m; if the wrong side is cut, flip the normal to 0 0 -1
        UImanager->ApplyCommand("/vis/viewer/zoomTo 2"); //scene radius is about 34 m, so about 17 m each side of the target
        UImanager->ApplyCommand("/vis/viewer/set/autorefresh true"); //changes for new run
        UImanager->ApplyCommand("/vis/scene/add/trajectories smooth"); //draws particle traj
        UImanager->ApplyCommand("/vis/scene/endOfEventAction accumulate"); //accumulates events to 1 run, not rlly needed
        UImanager->ApplyCommand("/vis/scene/add/scale 8 m z 1 1 1 manual 0 -2 33 m"); //white 8 m scale (bar, ticks and label take this colour) along z over the last 8 m (29 to 37 m, default zdcFrontZ), 2 m below the beam axis so it clears the 1 m wide blocks (to -0.5 m)
        UImanager->ApplyCommand("/vis/scene/add/scale 1 m y 1 1 1 manual 0 0 38.6 m"); //white 1 m vertical scale, the height of the world box (y -0.5 to 0.5 m), just outside the box end (box ends at z = 38.1 m with the default 37 m layout)
        //UImanager->ApplyCommand("/vis/scene/add/axes"); //axis for rotations @ center
        UImanager->ApplyCommand("/vis/scene/add/eventID"); //should which event happened
        UImanager->ApplyCommand("/run/printProgress 100"); //one progress line every 100 events
        UImanager->ApplyCommand("/tracking/verbose 0"); //step log off in the GUI, set /tracking/verbose 1 in the session to get it back
    
        //UImanager->ApplyCommand("/vis/filtering/trajectories/create/particleFilter");
        //UImanager->ApplyCommand("/vis/filtering/trajectories/particleFilter-0/add proton");
    } else{
        G4String command = "/control/execute ";
        G4String fileName = argv[1];
        UImanager->ApplyCommand(command+fileName);
    }; //when you use ./sim [file].mac, runs file.mac

    if (ui != nullptr){
    ui->SessionStart();
    delete ui;
    }

    delete runManager;
    return 0;
}
