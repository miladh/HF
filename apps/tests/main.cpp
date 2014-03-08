#include <unittest++/UnitTest++.h>
#include <unittest++/Test.h>
#include <unittest++/TestReporterStdout.h>
#include <unittest++/TestRunner.h>

#include <hf.h>

int main()
{
    int result = 0;
    bool runAllTests = 1;

    UnitTest::TestReporterStdout reporter;
    UnitTest::TestRunner runner(reporter);


    result += runner.RunTestsIf(UnitTest::Test::GetTestList(), "DEVELOPMENT", UnitTest::True(), 0);

    if(runAllTests){
        result += runner.RunTestsIf(UnitTest::Test::GetTestList(), "SLOWTESTS", UnitTest::True(), 0);
    }

    return result;
}
