#include "catch2/catch.hpp"

//#include "share/scream_types.hpp"
#include <algorithm>
#include <array>
#include <random>
#include <thread>

#include "ekat/scream_kokkos.hpp"
#include "ekat/scream_pack.hpp"
#include "ekat/scream_types.hpp"
#include "ekat/util/scream_arch.hpp"
#include "ekat/util/scream_kokkos_utils.hpp"
#include "ekat/util/scream_utils.hpp"
#include "physics/share/physics_constants.hpp"
#include "physics/shoc/shoc_functions.hpp"
#include "physics/shoc/shoc_functions_f90.hpp"
#include "shoc_unit_tests_common.hpp"

namespace scream {
namespace shoc {
namespace unit_test {

template <typename D>
struct UnitWrap::UnitTest<D>::TestShocPdfComputeQs {

  static void run_property()
  {
  
    // Property tests for the SHOC function
    //  shoc_assumed_pdf_compute_qs

    // TEST ONE
    // Zero Skewness tests.  Given inputs where Temperatures for
    //  the two gaussians are the same, verify that outputs are the same
    
    // Define temperature of Gaussian 1 [K]
    static constexpr Real Tl1_1_eq = 290.0;
    // Define temperature of Gaussian 2 [K]
    static constexpr Real Tl1_2_eq = Tl1_1_eq;
    // Define pressure value [Pa]
    static constexpr Real pval_eq = 85000.0;
    
    // Initialize data structure for bridging to F90
    SHOCPDFcompqsData SDS;
    
    // Fill in input data
    SDS.Tl1_1 = Tl1_1_eq;
    SDS.Tl1_2 = Tl1_2_eq;
    SDS.pval = pval_eq;
    
    // Check the inputs
    REQUIRE(SDS.Tl1_1 == SDS.Tl1_2);
    REQUIRE(SDS.pval > 0);
    
    // Call the fortran implementation 
    shoc_assumed_pdf_compute_qs(SDS);
    
    // Check the result
    REQUIRE(SDS.beta1 == SDS.beta2);
    REQUIRE(SDS.qs1 == SDS.qs2);
    
    // TEST TWO
    // Pressure test.  Use the data from the last test and modify
    //  the pressure level to be lower and verify that qs1 is lower.
    
    static constexpr Real pval_high = 100000.0;
    
    // Save the result from last test
    Real qs1_test1 = SDS.qs1;
    
    // Load new input
    SDS.pval = pval_high;
    
    // Verify new pressure is greater than last test
    REQUIRE(SDS.pval > pval_eq);
    
    // Call the fortran implementation 
    shoc_assumed_pdf_compute_qs(SDS);
    
    // Check the result
    REQUIRE(SDS.qs1 < qs1_test1);  
    
    // TEST THREE
    // Skewed test.  using pressre input from last test, but using
    //  a skewed temperature distribution.  Verify that the gaussian
    //  with lower input temperature also has lower qs value
    
    // Define temperature of Gaussian 1 [K]
    static constexpr Real Tl1_1_skew = 290;
    // Define temperature of Gaussian 2 [K]
    static constexpr Real Tl1_2_skew = 300;      
    
    // Load new input
    SDS.Tl1_1 = Tl1_1_skew;
    SDS.Tl1_2 = Tl1_2_skew;
    
    // Verify Gaussians are not equal
    REQUIRE(SDS.Tl1_1 != SDS.Tl1_2);
    
    // Call the fortran implementation 
    shoc_assumed_pdf_compute_qs(SDS);  
    
    // Check the result
    if (SDS.Tl1_1 < SDS.Tl1_2){
      REQUIRE(SDS.qs1 < SDS.qs2);
      REQUIRE(SDS.beta1 > SDS.beta2);
    }  
    else{
      REQUIRE(SDS.qs1 > SDS.qs2);
      REQUIRE(SDS.beta1 < SDS.beta2);
    }

  }
  
};

}  // namespace unit_test
}  // namespace shoc
}  // namespace scream

namespace {

TEST_CASE("shoc_pdf_compute_qs_property", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestShocPdfComputeQs;

  TestStruct::run_property();
}

TEST_CASE("shoc_pdf_compute_qs_b4b", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestShocPdfComputeQs;

}

} // namespace
