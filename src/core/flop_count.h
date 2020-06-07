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
    FlopCount () : flop_sum{0} { ; }
    static FlopCount without_instr_mix(double flops) { FlopCount fc; fc.flop_sum = flops; return fc; }
    double flop_sum;
    std::map<std::string, size_t> instr_mix;

    FlopCount& operator+=(const FlopCount& rhs) { 
        this->flop_sum += rhs.flop_sum;
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
        this->flop_sum += instr.flop_cost;
        this->instr_mix[instr.instr_name]++;
        return *this;
    };
    FlopCount operator+(const Instruction& rhs) {
        FlopCount lhs = *this;
        lhs += rhs;
        return lhs;
    }
};

inline FlopCount operator*(const size_t& lhs, const FlopCount& rhs) {
    FlopCount ret;
    for (size_t i = 0; i < lhs; i++) {
        ret += rhs;
    }
    return ret;
}

inline FlopCount operator*(const double& lhs, const Instruction& rhs) {
    FlopCount ret;
    ret += rhs;
    ret.flop_sum *= lhs;
    ret.instr_mix.begin()->second *= lhs;
    return ret;
}

inline FlopCount operator+(const Instruction& lhs, const Instruction& rhs) {
    FlopCount ret;
    ret += lhs;
    ret += rhs;
    return ret;
}

inline FlopCount operator+(const Instruction& lhs, FlopCount rhs) {
    FlopCount ret = rhs;
    ret += lhs;
    return ret;
}


class AbsInstr:        public Instruction { public: AbsInstr():        Instruction("abs", 0.0)  {}}; // TODO
class SinInstr:        public Instruction { public: SinInstr():        Instruction("sin", 15.0) {}}; //All values set to 1.0 for now, we change later
class Atan2Instr:      public Instruction { public: Atan2Instr():      Instruction("atan2", 25.0) {}};
class MulInstr:        public Instruction { public: MulInstr():        Instruction("mul", 1.0)  {}};
class AddInstr:        public Instruction { public: AddInstr():        Instruction("add", 1.0)  {}};
class DivInstr:        public Instruction { public: DivInstr():        Instruction("div", 5.0)  {}};
class SqrtInstr:       public Instruction { public: SqrtInstr():       Instruction("sqrt", 7.0)  {}};
class RsqrtInstr:      public Instruction { public: RsqrtInstr():      Instruction("rsqrt", 1.0)  {}};
class NegationInstr:   public Instruction { public: NegationInstr():   Instruction("negation", 1.0)  {}};
class CosInstr:        public Instruction { public: CosInstr():        Instruction("cos", 15.0) {}};
class DoublecompInstr: public Instruction { public: DoublecompInstr(): Instruction("doublecomp", 1.0)  {}};
class RandInstr:       public Instruction { public: RandInstr():       Instruction("rand", 30.0) {}}; //Guess
class FastrandInstr:   public Instruction { public: FastrandInstr():   Instruction("fastrand", 7.0)  {}}; //Guess
class FloorInstr:      public Instruction { public: FloorInstr():      Instruction("floor", 1.0)  {}};
class ModuloInstr:     public Instruction { public: ModuloInstr():     Instruction("modulo", 15.0) {}}; //Guess
class PowInstr:        public Instruction { public: PowInstr():        Instruction("pow", 1.0)  {}}; // 1 Mult  = pow2 (The only one we use)
class ExpInstr:        public Instruction { public: ExpInstr():        Instruction("exp", 10.0) {}};
