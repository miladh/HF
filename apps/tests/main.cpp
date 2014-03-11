#include <unittest++/UnitTest++.h>
#include <unittest++/Test.h>
#include <unittest++/TestReporterStdout.h>
#include <unittest++/TestRunner.h>
#include <mpi.h>
#include <hf.h>

int main(int argc, char **argv)
{
    MPI::Init(argc, argv);
    int result = 0;
    bool slowTests = 1;
    bool gradient = 1;

    UnitTest::TestReporterStdout reporter;
    UnitTest::TestRunner runner(reporter);

 if(MPI::COMM_WORLD.Get_rank()==0){
    result += runner.RunTestsIf(UnitTest::Test::GetTestList(), "DEVELOPMENT", UnitTest::True(), 0);

    if(slowTests){
        result += runner.RunTestsIf(UnitTest::Test::GetTestList(), "SLOWTESTS", UnitTest::True(), 0);
    }
    if(gradient){
        result += runner.RunTestsIf(UnitTest::Test::GetTestList(), "GRADIENT", UnitTest::True(), 0);

    }

    return result;
 }
}
