#include <iostream>
#include <vector>
#include <string>
#include <chrono>
#include <array>

#define N 10000000
unsigned int composite_rule_r = 3;

// create bit masks of length 2r + 1 for truncating configurations
constexpr std::array<uint64_t, 30>  create_masks() {
    std::array<uint64_t, 30> table = {0};
    for (unsigned int i = 1; i < 30; ++i) {
        table[i] = (1UL << (2*i + 1)) - 1;
    }
    return table;
}

constexpr std::array<uint64_t, 30> masks = create_masks();

std::vector<unsigned char> config;

unsigned int config_vec_size;
unsigned int config_width;
unsigned int i;
unsigned int tau;

class Timer {
private:
    using Clock = std::chrono::steady_clock;
    using Second = std::chrono::duration<double, std::ratio<1>>;
    std::chrono::time_point<Clock> m_beg { Clock::now() };

public:
    void reset() {
        m_beg = Clock::now();
    }

    double elapsed() const {
        return std::chrono::duration_cast<Second>(Clock::now() - m_beg).count();
    }
};

void pad_config(){
    config.insert(config.begin(), 6000, 0);
    config.insert(config.end(), 6000, 0);
    config_vec_size = config.size();
}

inline void run_simulation() {

    const unsigned int r = 1;

    Timer t;
    unsigned char config_ptr;

    for (; tau < N; ++tau) {
        
        unsigned int j;
        
        if (2*(config_width + 2*r + 1) >= config_vec_size)
            pad_config();

        config_ptr=0;
        for (j = config_vec_size / 2 - config_width; j < config_vec_size / 2 + config_width; j++){
            config_ptr <<= 1;
            config_ptr |= config[j];
            //explicit boolean equation for rule 30 in terms of XOR and OR
            config[j-1] = ((config_ptr>>2)&1)^(((config_ptr>>1)&1) | (config_ptr&1));
        }
        
        config_width += r;   
    
        if (tau % 500 == 0)
            std::cout << t.elapsed() << '\t' << tau << std::endl;
    }
}

inline void run_precomputed_simulation() {

    //rule 30 lookup
    std::vector<unsigned char> init_rule = {0,1,1,1,1,0,0,0};
    
    std::vector<unsigned char> composite_rule;

    Timer t;
    
    std::cout << t.elapsed() << '\t' << tau << std::endl;
    
    unsigned int j;
    unsigned int k;
    unsigned int l;

    //precompute the composite rule starting from r = 1
    const unsigned int r = composite_rule_r;
    uint64_t local_config[2] = {0,0};
    composite_rule.resize(1UL << (2*r + 1));

    for (j = 0; j < (1UL << (2*r+ 1)); ++j){ 
        
        local_config[0] = j;

        for (k  = 0; k < r; ++k){
            
            local_config[(k+1)%2] = 0;

            for (l = 0; l < 2*r + 1 - 2*k; ++l){
                local_config[(k+1)%2] |= static_cast<uint64_t>(init_rule[(local_config[k%2] >> l) & 0b111]) << l;
            }
        }

        composite_rule[j] = local_config[k%2] & 1UL;
    }

    //give sufficient padding when starting from a simple seed
    config_width = 3*r;

    std::cout << t.elapsed() << '\t' << tau << std::endl;
    uint64_t config_ptr = 0;

    //run the simulation
    for (; tau < N; ++i, tau += r){

        if (2*(config_width + (2*r + 1)) >= config_vec_size)
            pad_config();
        
        config_ptr = 0;
        for (j = config_vec_size/2 - config_width; j < config_vec_size/2 + config_width; j++) {
            config_ptr <<= 1;
            config_ptr |= config[j];
            config_ptr &= masks[r];
            config[j - r] = composite_rule[config_ptr];            
        }

        config_width += r;   

        if (i % (500/r) == 0)
            std::cout << t.elapsed() << '\t' << tau << std::endl;
    }
}

int main(int argc, char* argv[]) {

    if (argc < 2){
        std::cout << "Enter a composite rule radius or 0 for the reference: " << std::endl;
        std::string input;
        std::cin >> input;
        composite_rule_r = std::stoi(input);
    }
    else
        composite_rule_r = std::stoi(argv[1]);
        
    i = 1;
    tau = 1;
        
    config.resize(6000);
    config_vec_size = config.size();

    //initialize the config with a simple seed
    config[config_vec_size/2] = 1;

    if (composite_rule_r == 0)
        run_simulation();
    else
        run_precomputed_simulation();    

    return 0;

}
