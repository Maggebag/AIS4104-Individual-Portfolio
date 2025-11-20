#include <utility/enum.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

using namespace AIS4104;

enum class UnitTestEnum
{
    ALPHA, BRAVO, CHARLIE, DELTA
};

static std::array<UnitTestEnum, 4> enum_options = {
    UnitTestEnum::ALPHA,
    UnitTestEnum::BRAVO,
    UnitTestEnum::CHARLIE,
    UnitTestEnum::DELTA
};
static std::array<const char*, 4> enum_labels = {"ALPHA", "BRAVO", "CHARLIE", "DELTA"};

utility::Enum<UnitTestEnum, 4> unit_test_enum(UnitTestEnum initial_value = UnitTestEnum::DELTA)
{
    return utility::Enum<UnitTestEnum, 4>(
        initial_value,
        {
            std::make_pair(UnitTestEnum::ALPHA, enum_labels[0]),
            std::make_pair(UnitTestEnum::BRAVO, enum_labels[1]),
            std::make_pair(UnitTestEnum::CHARLIE, enum_labels[2]),
            std::make_pair(UnitTestEnum::DELTA, enum_labels[3])
        }
    );
}

TEST_CASE("enum_options_equals_in_value_size_and_order")
{
    auto options = unit_test_enum().options();
    REQUIRE(options.size() == enum_options.size());

    for(auto i = 0u; i < enum_options.size(); i++)
        REQUIRE(options[i] == enum_options[i]);
}

TEST_CASE("enum_labels_equals_in_value_size_and_order")
{
    auto labels = unit_test_enum().labels();
    REQUIRE(labels.size() == enum_labels.size());

    for(auto i = 0u; i < enum_labels.size(); i++)
        REQUIRE(std::string(labels[i]) == enum_labels[i]);
}

TEST_CASE("enum_values_and_labels_map_correctly")
{
    auto test_enum = unit_test_enum();

    for(auto i = 0u; i < enum_options.size(); i++)
        REQUIRE(test_enum.label(enum_options[i]) == enum_labels[i]);
}

TEST_CASE("enum_initial_value_and_label_set")
{
    auto test_idx = 2u;
    auto test_enum = unit_test_enum(enum_options[test_idx]);

    REQUIRE(test_enum.value() == enum_options[test_idx]);
    REQUIRE(test_enum.label() == std::string(enum_labels[test_idx]));
}

TEST_CASE("enum_assignment_operator_correctly_updates_value_and_label")
{
    auto initial_idx = 2u;
    auto assigned_idx = 3u;
    auto test_enum = unit_test_enum(enum_options[initial_idx]);

    REQUIRE(test_enum.value() == enum_options[initial_idx]);
    REQUIRE(test_enum.label() == std::string(enum_labels[initial_idx]));

    test_enum = enum_options[assigned_idx];
    REQUIRE(test_enum.value() == enum_options[assigned_idx]);
    REQUIRE(test_enum.label() == std::string(enum_labels[assigned_idx]));
}

TEST_CASE("enum_assignment_from_value_correctly_updates_value_and_label")
{
    auto initial_idx = 2u;
    auto assigned_idx = 3u;
    auto test_enum = unit_test_enum(enum_options[initial_idx]);

    REQUIRE(test_enum.value() == enum_options[initial_idx]);
    REQUIRE(test_enum.label() == std::string(enum_labels[initial_idx]));

    test_enum.set(enum_options[assigned_idx]);
    REQUIRE(test_enum.value() == enum_options[assigned_idx]);
    REQUIRE(test_enum.label() == std::string(enum_labels[assigned_idx]));
}

TEST_CASE("enum_assignment_from_label_correctly_updates_value_and_label")
{
    auto initial_idx = 2u;
    auto assigned_idx = 3u;
    auto test_enum = unit_test_enum(enum_options[initial_idx]);

    REQUIRE(test_enum.value() == enum_options[initial_idx]);
    REQUIRE(test_enum.label() == std::string(enum_labels[initial_idx]));

    test_enum.set(enum_labels[assigned_idx]);
    REQUIRE(test_enum.value() == enum_options[assigned_idx]);
    REQUIRE(test_enum.label() == std::string(enum_labels[assigned_idx]));
}

TEST_CASE("enum_equality_operators_for_values")
{
    auto initial_idx = 2u;
    auto test_enum = unit_test_enum(enum_options[initial_idx]);

    REQUIRE(test_enum == enum_options[initial_idx]);

    REQUIRE(test_enum != enum_options[0u]);
    REQUIRE(test_enum != enum_options[1u]);
    REQUIRE(test_enum != enum_options[3u]);
}

TEST_CASE("enum_equality_operators_for_labels")
{
    auto initial_idx = 2u;
    auto test_enum = unit_test_enum(enum_options[initial_idx]);

    REQUIRE(test_enum == enum_labels[initial_idx]);

    REQUIRE(test_enum != enum_labels[0u]);
    REQUIRE(test_enum != enum_labels[1u]);
    REQUIRE(test_enum != enum_labels[3u]);
}
