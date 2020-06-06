#pragma once
#include <string>
#include <map>

class Instruction {
public:
    Instruction(std::string name_, double cost_) : instr_name{name_}, flop_cost{cost_} {};
    std::string const instr_name;
    double const flop_cost;
};

class FlopCount {
public:
    double flop_count;
    std::map<std::string, size_t> instr_mix;
    FlopCount& operator+=(const FlopCount& rhs) { 
        this->flop_count += rhs.flop_count;
        for (auto key_val : rhs.instr_mix) {
            this->instr_mix[key_val.first] += key_val.second;
        }
        return *this;
    };
    FlopCount operator+(const FlopCount& rhs) {
        FlopCount lhs = *this;
        lhs += rhs;
        return lhs;
    }

    FlopCount& operator+=(const Instruction& instr) { 
        this->flop_count += instr.flop_cost;
        this->instr_mix[instr.instr_name]++;
        return *this;
    };
    FlopCount operator+(const Instruction& rhs) {
        FlopCount lhs = *this;
        lhs += rhs;
        return lhs;
    }
};

FlopCount operator*(const double& lhs, const Instruction& rhs) {
    FlopCount ret;
    ret += rhs;
    ret.flop_count *= lhs;
    ret.instr_mix.begin()->second *= lhs;
    return ret;
}

FlopCount operator+(const Instruction& lhs, const Instruction& rhs) {
    FlopCount ret;
    ret += lhs;
    ret += rhs;
    return ret;
}


class AbsInstr:        public Instruction { AbsInstr():        Instruction("abs", 0.0)  {}}; // TODO
class SinInstr:        public Instruction { SinInstr():        Instruction("sin", 15.0) {}}; //All values set to 1.0 for now, we change later
class Atan2Instr:      public Instruction { Atan2Instr():      Instruction("atan2", 25.0) {}};
class MulInstr:        public Instruction { MulInstr():        Instruction("mul", 1.0)  {}};
class AddInstr:        public Instruction { AddInstr():        Instruction("add", 1.0)  {}};
class DivInstr:        public Instruction { DivInstr():        Instruction("div", 5.0)  {}};
class SqrtInstr:       public Instruction { SqrtInstr():       Instruction("sqrt", 7.0)  {}};
class RsqrtInstr:      public Instruction { RsqrtInstr():      Instruction("rsqrt", 1.0)  {}};
class NegationInstr:   public Instruction { NegationInstr():   Instruction("negation", 1.0)  {}};
class CosInstr:        public Instruction { CosInstr():        Instruction("cos", 15.0) {}};
class DoublecompInstr: public Instruction { DoublecompInstr(): Instruction("doublecomp", 1.0)  {}};
class RandInstr:       public Instruction { RandInstr():       Instruction("rand", 30.0) {}}; //Guess
class FastrandInstr:   public Instruction { FastrandInstr():   Instruction("fastrand", 7.0)  {}}; //Guess
class FloorInstr:      public Instruction { FloorInstr():      Instruction("floor", 1.0)  {}};
class ModuloInstr:     public Instruction { ModuloInstr():     Instruction("modulo", 15.0) {}}; //Guess
class PowInstr:        public Instruction { PowInstr():        Instruction("pow", 1.0)  {}}; // 1 Mult  = pow2 (The only one we use)
class ExpInstr:        public Instruction { ExpInstr():        Instruction("exp", 10.0) {}};
