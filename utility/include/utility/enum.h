#ifndef UTILITY_ENUM_H
#define UTILITY_ENUM_H

#include <array>
#include <string>
#include <stdint.h>
#include <stdexcept>

namespace AIS4104::utility {

template<typename E, uint8_t N>
class Enum
{
public:
    Enum(E initial_value, const std::array<std::pair<E, const char*>, N> &labels)
        : m_value(initial_value)
    {
        for(auto i = 0u; i < N; i++)
        {
            m_values[i] = labels[i].first;
            m_labels[i] = labels[i].second;
        }
    }

    uint8_t value_index() const
    {
        for(auto i = 0u; i < N; i++)
            if(m_value == m_values[i])
                return i;
        return 0u;
    }

    std::array<E, N> options() const
    {
        return m_values;
    }

    const std::array<const char*, N>& labels() const
    {
        return m_labels;
    }

    E value() const
    {
        return m_value;
    }

    E value(const std::string &label)
    {
        for(auto i = 0u; i < N; i++)
            if(label == m_labels[i])
                return m_values[i];
        throw std::invalid_argument("Invalid enum Label");
    }

    std::string label() const
    {
        for(auto i = 0u; i < N; i++)
            if(m_values[i] == m_value)
                return m_labels[i];
        return std::string();
    }

    std::string label(E value) const
    {
        for(auto i = 0u; i < N; i++)
            if(m_values[i] == value)
                return m_labels[i];
        return std::string();
    }

    void set(E rhs)
    {
        m_value = rhs;
    }

    void set(const std::string &label)
    {
        m_value = value(label);
    }

    void operator=(E rhs)
    {
        m_value = rhs;
    }

    bool operator==(E rhs)
    {
        return m_value == rhs;
    }

    bool operator!=(E rhs)
    {
        return m_value != rhs;
    }

    bool operator==(const std::string &rhs)
    {
        return m_value == value(rhs);
    }

    bool operator!=(const std::string &rhs)
    {
        return m_value != value(rhs);
    }

private:
    E m_value;
    std::array<E, N> m_values;
    std::array<const char*, N> m_labels;
};

}

#endif
