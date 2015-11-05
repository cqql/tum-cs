#include <iostream>
#include <string>

class Person {
public:
  std::string name;
  std::string country;

  virtual void setName(std::string name) = 0;
  virtual std::string getName() = 0;

  virtual void setCountry(std::string country) = 0;
  virtual std::string getCountry() = 0;
private:
  int age;
};

class Student : public Person {
public:
  Student(std::string name, std::string country) {
    this->name = name;
    this->country = country;
  }

  void setName(std::string name) {
    this->name = name;
  }

  std::string getName() {
    return this->name;
  }

  void setCountry(std::string country) {
    this->country = country;
  }

  std::string getCountry() {
    return this->country;
  }
};

int main() {
  Student student("Marten", "de");

  std::cout << student.getName() << std::endl;

  return 0;
}
