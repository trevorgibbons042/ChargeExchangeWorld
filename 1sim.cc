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
#include "G4HadronicParameters.hh"

int main(int argc, char** argv)
{ 
    G4RunManager *runManager = new G4RunManager(); //initalizes for all
    runManager->SetUserInitialization(new MyDetectorConstruction()); //construction.cc stuff
    runManager->SetUserInitialization(new MyPhysicsList()); //physics.cc stuff
    runManager->SetUserInitialization(new MyActionInitialization()); //action.cc stuff
    //no need to load generator as its with action

    runManager->Initialize(); //STARTS IT UP!, get rid of particles arent there
    
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
        UImanager->ApplyCommand("/vis/viewer/set/autorefresh true"); //changes for new run
        UImanager->ApplyCommand("/vis/scene/add/trajectories smooth"); //draws particle traj
        UImanager->ApplyCommand("/vis/scene/endOfEventAction accumulate"); //accumulates events to 1 run, not rlly needed
        UImanager->ApplyCommand("/vis/scene/add/scale 10 cm"); //adds a scale
        //UImanager->ApplyCommand("/vis/scene/add/axes"); //axis for rotations @ center
        UImanager->ApplyCommand("/vis/scene/add/eventID"); //should which event happened
        UImanager->ApplyCommand("/tracking/verbose 1");
    
        //UImanager->ApplyCommand("/vis/filtering/trajectories/create/particleFilter");
        //UImanager->ApplyCommand("/vis/filtering/trajectories/particleFilter-0/add proton");
    } else{
        G4String command = "/control/execute ";
        G4String fileName = argv[1];
        UImanager->ApplyCommand(command+fileName);
    }; //when you use ./sim [file].mac, runs file.mac

    if (ui != nullptr) {
    ui->SessionStart();
    delete ui;
    }
    delete runManager;
    return 0;
}
