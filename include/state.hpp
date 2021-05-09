#ifndef STATE_HPP
#define STATE_HPP

#include "defines.hpp"



class State {

    // ----------------- TRANSIENT STATES -----------------


    // ----------------- PERSISTENT STATES ----------------



    // ----------------------------------------------------

    static State* state;
    State(){}
public:

    [[nodiscard]] float deltaTime() const;
    [[nodiscard]] size_t fps() const;
    
};

class UpdateState : public State {
    public:

};

#endif //STATE_HPP