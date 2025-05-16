#!/bin/bash

# Colors for output
GREEN='\033[0;32m'
RED='\033[0;31m'
NC='\033[0m' # No Color

# Clean and build the project
echo "Building the project..."
make clean
make all

if [ $? -ne 0 ]; then
    echo -e "${RED}Build failed. Aborting tests.${NC}"
    exit 1
fi

echo -e "${GREEN}Build successful.${NC}"

# Directory for test results
mkdir -p test_results

# Run the main test suite
echo -e "\nRunning main test suite..."
./test > test_results/test_suite_output.txt
echo -e "${GREEN}Main test suite completed. Results in test_results/test_suite_output.txt${NC}"

# Function to run a specific test case and check results
run_test() {
    local a=$1
    local b=$2
    local c=$3
    local d=$4
    local nx=$5
    local ny=$6
    local k=$7
    local eps=$8
    local max_iter=$9
    local p=${10}
    local test_name=${11}

    echo "Running test: $test_name"
    ./a.out $a $b $c $d $nx $ny $k $eps $max_iter $p > test_results/${test_name}.txt

    if [ $? -eq 0 ]; then
        echo -e "${GREEN}Test $test_name completed successfully.${NC}"
    else
        echo -e "${RED}Test $test_name failed.${NC}"
    fi
}

# Test exactness for polynomials of degree 0 and 1
echo -e "\nRunning exactness tests..."
run_test 0 1 0 1 5 5 0 1e-10 1000 4 "exact_constant"
run_test 0 1 0 1 5 5 1 1e-10 1000 4 "exact_linear_x"
run_test 0 1 0 1 5 5 2 1e-10 1000 4 "exact_linear_y"
run_test 0 1 0 1 5 5 3 1e-10 1000 4 "exact_linear_xy"

# Test convergence rates
echo -e "\nRunning convergence tests..."
run_test 0 1 0 1 10 10 5 1e-8 1000 4 "convergence_10x10"
run_test 0 1 0 1 20 20 5 1e-8 1000 4 "convergence_20x20"
run_test 0 1 0 1 40 40 5 1e-8 1000 4 "convergence_40x40"

# Test all functions on a medium grid
echo -e "\nRunning function tests..."
for k in {0..7}; do
    run_test 0 1 0 1 20 20 $k 1e-6 1000 4 "function_k${k}"
done

# Test parallel scaling
echo -e "\nRunning parallel scaling tests..."
for p in 1 2 4 8; do
    run_test 0 1 0 1 30 30 6 1e-6 1000 $p "parallel_p${p}"
done

# Test different domain shapes
echo -e "\nRunning domain shape tests..."
run_test 0 1 0 1 20 20 5 1e-6 1000 4 "domain_square"
run_test 0 2 0 1 20 20 5 1e-6 1000 4 "domain_rectangle"
run_test -1 1 -1 1 20 20 5 1e-6 1000 4 "domain_negative"

# Test with different solver tolerances
echo -e "\nRunning solver tolerance tests..."
run_test 0 1 0 1 20 20 5 1e-4 1000 4 "tolerance_1e-4"
run_test 0 1 0 1 20 20 5 1e-6 1000 4 "tolerance_1e-6"
run_test 0 1 0 1 20 20 5 1e-8 1000 4 "tolerance_1e-8"

# Test with different maximum iterations
echo -e "\nRunning iteration limit tests..."
run_test 0 1 0 1 20 20 6 1e-6 100 4 "max_iter_100"
run_test 0 1 0 1 20 20 6 1e-6 500 4 "max_iter_500"
run_test 0 1 0 1 20 20 6 1e-6 1000 4 "max_iter_1000"

# Check polynomial exactness from results
echo -e "\nVerifying polynomial exactness..."
for test in exact_constant exact_linear_x exact_linear_y exact_linear_xy; do
    # Extract R3 error (max error at grid points)
    r3=$(grep "R3 =" test_results/${test}.txt | awk '{print $5}')

    # Check if error is small enough to be considered "exact"
    if (( $(echo "$r3 < 1e-8" | bc -l) )); then
        echo -e "${GREEN}Test $test PASSED: Error R3 = $r3 < 1e-8${NC}"
    else
        echo -e "${RED}Test $test FAILED: Error R3 = $r3 >= 1e-8${NC}"
    fi
done

# Check convergence rates
echo -e "\nVerifying convergence rates..."
r3_10=$(grep "R3 =" test_results/convergence_10x10.txt | awk '{print $5}')
r3_20=$(grep "R3 =" test_results/convergence_20x20.txt | awk '{print $5}')
r3_40=$(grep "R3 =" test_results/convergence_40x40.txt | awk '{print $5}')

ratio1=$(echo "scale=2; $r3_10 / $r3_20" | bc)
ratio2=$(echo "scale=2; $r3_20 / $r3_40" | bc)

echo "Error ratio 10x10 to 20x20: $ratio1 (expect ~4)"
echo "Error ratio 20x20 to 40x40: $ratio2 (expect ~4)"

if (( $(echo "$ratio1 > 3.0" | bc -l) )) && (( $(echo "$ratio2 > 3.0" | bc -l) )); then
    echo -e "${GREEN}Convergence test PASSED: Ratios are close to expected quadratic convergence${NC}"
else
    echo -e "${RED}Convergence test WARNING: Convergence rates may not be optimal${NC}"
fi

# Check parallel scaling
echo -e "\nVerifying parallel scaling..."
t1_p1=$(grep "T1 =" test_results/parallel_p1.txt | awk '{print $5}')
t1_p2=$(grep "T1 =" test_results/parallel_p2.txt | awk '{print $5}')
t1_p4=$(grep "T1 =" test_results/parallel_p4.txt | awk '{print $5}')
t1_p8=$(grep "T1 =" test_results/parallel_p8.txt | awk '{print $5}')

speedup_p2=$(echo "scale=2; $t1_p1 / $t1_p2" | bc)
speedup_p4=$(echo "scale=2; $t1_p1 / $t1_p4" | bc)
speedup_p8=$(echo "scale=2; $t1_p1 / $t1_p8" | bc)

echo "Speedup with 2 threads: $speedup_p2 (ideal: 2.0)"
echo "Speedup with 4 threads: $speedup_p4 (ideal: 4.0)"
echo "Speedup with 8 threads: $speedup_p8 (ideal: 8.0)"

# Final summary
echo -e "\n${GREEN}All tests completed. Check test_results/ directory for detailed outputs.${NC}"
