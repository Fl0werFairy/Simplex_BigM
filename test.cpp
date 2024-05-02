#include <iostream>
#include <gtest/gtest.h>
#include "modified.h"
#include "golden.h"
#include <fstream>
#include <array>

using namespace std;
using namespace std::literals;
array input_files = {
		"A_unbounded.txt"s,
		"A2.txt"s,
		"B1.txt"s,
		"B_infeasible.txt"s
};

class enumTest : public testing::TestWithParam<string>
{
};

TEST_P(enumTest, equivelantness)
{
	string input_path = GetParam();
	cout << "file_name: " << input_path << '\n';
	ifstream fin1{input_path}, fin2{input_path};
	auto [kg, ansg] = golden::work(fin1);
	auto [kt, anst] = trial::work(fin2);
	EXPECT_DOUBLE_EQ(kg, kt);
	ASSERT_EQ(ansg.size(), anst.size());
	for (int i = 0; i < ansg.size(); ++i)EXPECT_DOUBLE_EQ(ansg[i], anst[i]);
}

INSTANTIATE_TEST_SUITE_P(MyTests, enumTest, ::testing::ValuesIn(input_files));

int main()
{
	testing::InitGoogleTest();
	return RUN_ALL_TESTS();
}

