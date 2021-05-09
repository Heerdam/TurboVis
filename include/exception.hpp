
#include <exception>
#include <string>

class TVexception : public std::exception {
    const std::string msg;
 public:
    TVexception(const std::string& _msg) : msg(_msg){}
    virtual const char* what() const noexcept override {
        return msg.c_str();
    }
};