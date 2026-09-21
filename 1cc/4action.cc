#include "4action.hh"

MyActionInitialization::MyActionInitialization()
{};

MyActionInitialization::~MyActionInitialization()
{};

void MyActionInitialization::Build() const
{
    MyPrimaryGenerator *generator = new MyPrimaryGenerator(); //setup for 5generator.hh
    SetUserAction(generator);

    MyRunAction *runAction = new MyRunAction(); //setup for 7run.hh
    SetUserAction(runAction);

    MyEventAction *eventAction = new MyEventAction(runAction); //setup for 8event.hh
    SetUserAction(eventAction);

    MySteppingAction *steppingAction = new MySteppingAction(eventAction); //setup for 8event.hh (part 2)
    SetUserAction(steppingAction);
};