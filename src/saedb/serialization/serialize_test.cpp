#include "serialization_includes.hpp"
#include "../test/testharness.hpp"
#include <fstream>
#include <cstring>
#include <cstdio>
#include <list>
#include <vector>

using namespace std;
using namespace sae::serialization;

struct SerializationTest {
    SerializationTest() {

    }

    ~SerializationTest() {
        // clean all bin file
        remove("int.bin");
        remove("char-star.bin");
        remove("string.bin");
        remove("char-array.bin");
        remove("list.bin");
        remove("vector.bin");
    }
};

TEST(SerializationTest, SaveAndLoad) {
    int i = 42;
    int j;
    {
        std::ofstream fout("int.bin", std::fstream::binary);
        OSerializeStream encoder(&fout);
        encoder << i;
        fout.close();
    }

    {
        std::ifstream fin("int.bin", std::fstream::binary);
        ISerializeStream decoder(&fin);
        decoder >> j;
        ASSERT_EQ(i, j);
        fin.close();
    }
}

TEST(SerializationTest, BasicTypes) {
    const char* a = "hello world";
    const char* b;
    {
        std::ofstream fout("char-star.bin", std::fstream::binary);
        OSerializeStream encoder(&fout);
        encoder << a;
        fout.close();
    }

    {
        std::ifstream fin("char-star.bin", std::fstream::binary);
        ISerializeStream decoder(&fin);
        decoder >> b;
        ASSERT_TRUE(memcmp(a, b, sizeof(a)) == 0);
    }

    string s("source");
    string t;

    {
        std::ofstream fout("string.bin", std::fstream::binary);
        OSerializeStream encoder(&fout);
        encoder << s;
    }

    {
        std::ifstream fin("string.bin", std::fstream::binary);
        ISerializeStream decoder(&fin);
        decoder >> t;
        ASSERT_TRUE(memcmp(s.c_str(), t.c_str(), sizeof(s.c_str())) == 0);
    }

    char c[] = "hello world";
    char d[0];
    {
        std::ofstream fout("char-array.bin", std::fstream::binary);
        OSerializeStream encoder(&fout);
        encoder << c;
        fout.close();
    }

    {
        std::ifstream fin("char-array.bin", std::fstream::binary);
        ISerializeStream decoder(&fin);
        decoder >> d;
        ASSERT_TRUE(memcmp(c, d, sizeof(c)) == 0);
    }
}

TEST(SerializationTest, List) {
    list<int> s;
    list<int> t;
    s.push_back(1);
    s.push_back(2);

    {
        ofstream fout("list.bin", std::fstream::binary);
        OSerializeStream encoder(&fout);
        encoder << s;
        fout.close();
    }

    {
        ifstream fin("list.bin", fstream::binary);
        ISerializeStream decoder(&fin);
        decoder >> t;
        ASSERT_EQ(t.size(), s.size());
        auto sbegin = s.begin();
        auto tbegin = t.begin();
        for (size_t i = 0;i < s.size(); ++i) {
            ASSERT_EQ(*sbegin, *tbegin);
            sbegin++;
            tbegin++;
        }
    }
}

TEST(SerializationTest, Vector) {
    vector<int> s;
    vector<int> t;
    s.push_back(1);
    s.push_back(2);

    {
        ofstream fout("vector.bin", fstream::binary);
        OSerializeStream encoder(&fout);
        encoder << s;
        fout.close();
    }

    {
        ifstream fin("vector.bin", fstream::binary);
        ISerializeStream decoder(&fin);
        decoder >> t;
        ASSERT_EQ(t.size(), s.size());
        auto sbegin = s.begin();
        auto tbegin = t.begin();
        for (size_t i = 0;i < s.size(); ++i) {
            ASSERT_EQ(*sbegin, *tbegin);
            sbegin++;
            tbegin++;
        }
    }
}


int main(){
    return ::saedb::test::RunAllTests();
}
