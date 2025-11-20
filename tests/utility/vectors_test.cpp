#include <utility/vectors.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "utility/math.h"

using namespace AIS4104;

Eigen::VectorXf test_vectorf()
{
    Eigen::VectorXf vectorf(4);
    vectorf << 2.f, 1.f, 6.f, 9.f;
    return vectorf;
}

Eigen::VectorXd test_vectord()
{
    Eigen::VectorXd vectord(4);
    vectord << 2.0, 1.0, 6.0, 9.0;
    return vectord;
}

std::vector<float> test_std_vectorf()
{
    return {2.f, 1.f, 6.f, 9.f};
}

std::vector<double> test_std_vectord()
{
    return {2.0, 1.0, 6.0, 9.0};
}

template<typename F> requires utility::is_floating_point<F>
F value_plus_eps(F a)
{
    return a + std::numeric_limits<F>::epsilon();
}

template<typename F> requires utility::is_floating_point<F>
F eps_factor(F n)
{
    return n * std::numeric_limits<F>::epsilon();
}

TEST_CASE("is_approx_equal_vectorf_std_vectorf")
{
    REQUIRE(utility::is_approx_equal(test_vectorf(), test_std_vectorf()));
    REQUIRE(!utility::is_approx_equal((value_plus_eps(1.f) * test_vectorf()).eval(), test_std_vectorf()));
}

TEST_CASE("is_approx_equal_vectord_std_vectorf")
{
    REQUIRE(utility::is_approx_equal(test_vectord(), test_std_vectorf()));
    REQUIRE(!utility::is_approx_equal((test_vectord() * value_plus_eps(1.f)).eval(), test_std_vectorf()));
}

TEST_CASE("is_approx_equal_vectord_std_vectord")
{
    REQUIRE(utility::is_approx_equal(test_vectord(), test_std_vectord()));
    REQUIRE(!utility::is_approx_equal((value_plus_eps(1.0) * test_vectord()).eval(), test_std_vectord()));
}

TEST_CASE("is_approx_equal_std_vectorf_vectorf")
{
    REQUIRE(utility::is_approx_equal(test_std_vectorf(), test_vectorf()));
    REQUIRE(!utility::is_approx_equal(test_std_vectorf(), (value_plus_eps(1.f)* test_vectorf()).eval()));
}

TEST_CASE("is_approx_equal_std_vectord_vectorf")
{
    REQUIRE(utility::is_approx_equal(test_std_vectorf(), test_vectord()));
    REQUIRE(!utility::is_approx_equal(test_std_vectorf(), (value_plus_eps(1.f) * test_vectord()).eval()));
}

TEST_CASE("is_approx_equal_std_vectord_vectord")
{
    REQUIRE(utility::is_approx_equal(test_std_vectord(), test_vectord()));
    REQUIRE(!utility::is_approx_equal(test_std_vectord(), (value_plus_eps(1.0) * test_vectord()).eval()));
}

TEST_CASE("is_approx_equal_kdl_vector_vectord")
{
    Eigen::Vector3d eigen_vector(2.0, 3.0, -1.0);

    KDL::Vector kdl_vector(2.0, 3.0, -1.0);

    REQUIRE(utility::is_approx_equal(kdl_vector, eigen_vector));
    REQUIRE(utility::is_approx_equal(eigen_vector, kdl_vector));

    REQUIRE(!utility::is_approx_equal(kdl_vector, (value_plus_eps(1.0) * eigen_vector).eval()));
    REQUIRE(!utility::is_approx_equal((value_plus_eps(1.0) * eigen_vector).eval(), kdl_vector));
}

TEST_CASE("is_approx_equal_kdl_rotation_matrix3d")
{
    Eigen::Vector3d zyx = Eigen::Vector3d{90.0, 45.0, 22.5} * utility::deg_to_rad;

    Eigen::Matrix3d eigen_rotation_1 = utility::rotation_matrix_from_euler_zyx(zyx);
    Eigen::Matrix3d eigen_rotation_2 = utility::rotation_matrix_from_euler_zyx((value_plus_eps(1.0) * zyx).eval());

    KDL::Rotation kdl_rotation = KDL::Rotation::EulerZYX(zyx.x(), zyx.y(), zyx.z());

    REQUIRE(utility::is_approx_equal(kdl_rotation, eigen_rotation_1));
    REQUIRE(utility::is_approx_equal(eigen_rotation_1, kdl_rotation));

    REQUIRE(!utility::is_approx_equal(kdl_rotation, eigen_rotation_2));
    REQUIRE(!utility::is_approx_equal(eigen_rotation_2, kdl_rotation));
}

TEST_CASE("is_approx_equal_kdl_frame_matrix4d")
{
    Eigen::Vector3d zyx = Eigen::Vector3d{90.0, 45.0, 22.5} * utility::deg_to_rad;

    Eigen::Vector3d eigen_vector(2.0, 3.0, -1.0);
    Eigen::Matrix3d eigen_rotation = utility::rotation_matrix_from_euler_zyx(zyx);
    Eigen::Matrix4d eigen_frame_1 = utility::transformation_matrix(eigen_rotation, eigen_vector);
    Eigen::Matrix4d eigen_frame_2 = utility::transformation_matrix(eigen_rotation, value_plus_eps(1.0) * eigen_vector);

    KDL::Vector kdl_vector(2.0, 3.0, -1.0);
    KDL::Rotation kdl_rotation = KDL::Rotation::EulerZYX(zyx.x(), zyx.y(), zyx.z());

    KDL::Frame kdl_frame = KDL::Frame(kdl_rotation, kdl_vector);

    REQUIRE(utility::is_approx_equal(kdl_frame, eigen_frame_1));
    REQUIRE(utility::is_approx_equal(eigen_frame_1, kdl_frame));

    REQUIRE(!utility::is_approx_equal(kdl_frame, eigen_frame_2));
    REQUIRE(!utility::is_approx_equal(eigen_frame_2, kdl_frame));
}

TEST_CASE("eigen_vectorf_to_std_vectorf")
{
    REQUIRE(utility::is_approx_equal(utility::to_std_vectorf(test_vectorf()), test_vectorf()));
    REQUIRE(!utility::is_approx_equal(utility::to_std_vectorf(test_vectorf()), (value_plus_eps(1.f) * test_vectorf()).eval()));
}

TEST_CASE("eigen_vectord_to_std_vectorf")
{
    REQUIRE(utility::is_approx_equal(utility::to_std_vectorf(test_vectord()), test_vectord()));
    REQUIRE(!utility::is_approx_equal(utility::to_std_vectorf(test_vectord()), (value_plus_eps(1.f) * test_vectord()).eval()));
}

TEST_CASE("eigen_vectord_to_std_vectord")
{
    REQUIRE(utility::is_approx_equal(utility::to_std_vectord(test_vectord()), test_vectord()));
    REQUIRE(!utility::is_approx_equal(utility::to_std_vectord(test_vectord()), (value_plus_eps(1.0) * test_vectord()).eval()));
}

TEST_CASE("eigen_vectorf_to_std_vectord")
{
    REQUIRE(utility::is_approx_equal(utility::to_std_vectord(test_vectorf()), test_vectorf()));
    REQUIRE(!utility::is_approx_equal(utility::to_std_vectord(test_vectorf()), (value_plus_eps(1.f) * test_vectorf()).eval()));
}

TEST_CASE("std_vectorf_to_eigen_vectord")
{
    REQUIRE(utility::is_approx_equal(utility::to_eigen_vectord(test_std_vectorf()), test_std_vectorf()));
    REQUIRE(!utility::is_approx_equal((value_plus_eps(1.0) * utility::to_eigen_vectord(test_std_vectorf())).eval(), test_std_vectorf()));
}

TEST_CASE("std_vectord_to_eigen_vectord")
{
    REQUIRE(utility::is_approx_equal(utility::to_eigen_vectord(test_std_vectord()), test_std_vectord()));
    REQUIRE(!utility::is_approx_equal((value_plus_eps(1.0) * utility::to_eigen_vectord(test_std_vectord())).eval(), test_std_vectord()));
}

TEST_CASE("std_vectorf_to_eigen_vectorf")
{
    REQUIRE(utility::is_approx_equal(utility::to_eigen_vectorf(test_std_vectorf()), test_std_vectorf()));
    REQUIRE(!utility::is_approx_equal((value_plus_eps(1.f) * utility::to_eigen_vectorf(test_std_vectorf())).eval(), test_std_vectorf()));
}

TEST_CASE("std_vectord_to_eigen_vectorf")
{
    REQUIRE(utility::is_approx_equal(utility::to_eigen_vectorf(test_std_vectord()), test_std_vectord()));
    REQUIRE(!utility::is_approx_equal((value_plus_eps(1.f) * utility::to_eigen_vectorf(test_std_vectord())).eval(), test_std_vectord()));
}

TEST_CASE("eigen_vectorf_to_eigen_vectord")
{
    REQUIRE(utility::is_approx_equal(utility::to_eigen_vectord(test_vectorf()), test_vectord()));
    REQUIRE(!utility::is_approx_equal(utility::to_eigen_vectord(test_vectorf()), (value_plus_eps(1.0) * test_vectord()).eval()));
}

TEST_CASE("eigen_vectord_to_eigen_vectorf")
{
    REQUIRE(utility::is_approx_equal(utility::to_eigen_vectorf(test_vectord()), test_vectorf()));
    REQUIRE(!utility::is_approx_equal(utility::to_eigen_vectorf(test_vectord()), (value_plus_eps(1.f) * test_vectorf()).eval()));
}

TEST_CASE("eigen_vectord_to_eigen_vectord")
{
    REQUIRE(utility::is_approx_equal(test_vectord(), test_vectord()));
    REQUIRE(!utility::is_approx_equal(test_vectord(), (value_plus_eps(1.0) * test_vectord()).eval()));
}

TEST_CASE("eigen_vectorf_to_eigen_vectorf")
{
    REQUIRE(utility::is_approx_equal(test_vectorf(), test_vectorf()));
    REQUIRE(!utility::is_approx_equal(test_vectorf(), (value_plus_eps(1.f) * test_vectorf()).eval()));
}

TEST_CASE("eigen_vectord_to_kdl_array")
{
    REQUIRE(utility::is_approx_equal(test_vectord(), utility::to_kdl_array(test_vectord()).data));
    REQUIRE(!utility::is_approx_equal(test_vectord(), utility::to_kdl_array(value_plus_eps(1.0) * test_vectord()).data));
}

TEST_CASE("kdl_vector_to_eigen_vector3d")
{
    KDL::Vector kdl_vector(2.0, 3.0, -1.0);
    Eigen::Vector3d eigen_vector(2.0, 3.0, -1.0);
    REQUIRE(utility::is_approx_equal(eigen_vector, utility::to_eigen_vector3d(kdl_vector)));
    REQUIRE(!utility::is_approx_equal(eigen_vector, (value_plus_eps(1.0) * utility::to_eigen_vector3d(kdl_vector)).eval()));
}

TEST_CASE("kdl_rotation_to_eigen_matrix3d")
{
    Eigen::Vector3d zyx = Eigen::Vector3d{90.0, 45.0, 22.5} * utility::deg_to_rad;
    Eigen::Matrix3d r = utility::rotation_matrix_from_euler_zyx(zyx);
    KDL::Rotation kdl_rotation = KDL::Rotation::EulerZYX(zyx.x(), zyx.y(), zyx.z());

    REQUIRE(utility::is_approx_equal(r, utility::to_eigen_matrix3d(kdl_rotation)));
}

TEST_CASE("kdl_frame_to_eigen_matrix4d")
{
    Eigen::Vector3d zyx = Eigen::Vector3d{90.0, 45.0, 22.5} * utility::deg_to_rad;

    Eigen::Vector3d eigen_vector(2.0, 3.0, -1.0);
    Eigen::Matrix3d eigen_rotation = utility::rotation_matrix_from_euler_zyx(zyx);
    Eigen::Matrix4d eigen_frame = utility::transformation_matrix(eigen_rotation, eigen_vector);

    KDL::Vector kdl_vector(2.0, 3.0, -1.0);
    KDL::Rotation kdl_rotation = KDL::Rotation::EulerZYX(zyx.x(), zyx.y(), zyx.z());

    KDL::Frame kdl_frame = KDL::Frame(kdl_rotation, kdl_vector);

    REQUIRE(utility::is_approx_equal(eigen_frame, utility::to_eigen_matrix4d(kdl_frame)));
}

TEST_CASE("eigen_vector3d_to_kdl_vector")
{
    KDL::Vector kdl_vector(2.0, 3.0, -1.0);
    Eigen::Vector3d eigen_vector(2.0, 3.0, -1.0);
    REQUIRE(utility::is_approx_equal(kdl_vector, utility::to_kdl_vector(eigen_vector)));
    REQUIRE(!utility::is_approx_equal(kdl_vector, utility::to_kdl_vector((value_plus_eps(1.0) *eigen_vector).eval())));
}

TEST_CASE("eigen_matrix3d_to_kdl_rotation")
{
    Eigen::Vector3d zyx = Eigen::Vector3d{90.0, 45.0, 22.5} * utility::deg_to_rad;
    Eigen::Matrix3d r = utility::rotation_matrix_from_euler_zyx(zyx);
    KDL::Rotation kdl_rotation = KDL::Rotation::EulerZYX(zyx.x(), zyx.y(), zyx.z());

    REQUIRE(utility::is_approx_equal(kdl_rotation, utility::to_kdl_rotation(r)));
}

TEST_CASE("eigen_matrix4d_to_kdl_frame")
{
    Eigen::Vector3d zyx = Eigen::Vector3d{90.0, 45.0, 22.5} * utility::deg_to_rad;

    Eigen::Vector3d eigen_vector(2.0, 3.0, -1.0);
    Eigen::Matrix3d eigen_rotation = utility::rotation_matrix_from_euler_zyx(zyx);
    Eigen::Matrix4d eigen_frame = utility::transformation_matrix(eigen_rotation, eigen_vector);

    KDL::Vector kdl_vector(2.0, 3.0, -1.0);
    KDL::Rotation kdl_rotation = KDL::Rotation::EulerZYX(zyx.x(), zyx.y(), zyx.z());

    KDL::Frame kdl_frame = KDL::Frame(kdl_rotation, kdl_vector);

    REQUIRE(utility::is_approx_equal(kdl_frame, utility::to_kdl_frame(eigen_frame)));
}
