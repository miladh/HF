#include <unittest++/UnitTest++.h>
#include <unittest++/Test.h>
#include <unittest++/TestReporterStdout.h>
#include <unittest++/TestRunner.h>
#include <boost/mpi.hpp>
#include <hf.h>

int main(int argc, char **argv)
{

#if USE_MPI
    boost::mpi::environment env(argc, argv);
#endif
    int result = 0;
    bool slowTests = 1;
    bool slowTests_UHF = 1;
    bool gradient = 1;

    UnitTest::TestReporterStdout reporter;
    UnitTest::TestRunner runner(reporter);

    result += runner.RunTestsIf(UnitTest::Test::GetTestList(), "DEVELOPMENT", UnitTest::True(), 0);

    if(slowTests){
        result += runner.RunTestsIf(UnitTest::Test::GetTestList(), "SLOWTESTS", UnitTest::True(), 0);
    }

    if(slowTests_UHF){
        result += runner.RunTestsIf(UnitTest::Test::GetTestList(), "SLOWTESTS_UHF", UnitTest::True(), 0);
    }

    if(gradient){
        result += runner.RunTestsIf(UnitTest::Test::GetTestList(), "GRADIENT", UnitTest::True(), 0);

    }

 return result;
}
