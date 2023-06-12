#include <iostream>
#include <list>

using namespace std;

struct Element {
    int m_val = -1;
    void set_val(int val);
    Element(int val): m_val(val) {}
};


int main() {
    /// @brief this code demonstrates using STL containers and their iterators

    // std::list implements a doubly-linked list
    list<Element> a = {Element(1), Element(2), Element(3), Element(5), Element(8)};

    // the iterator type associated to list<Element>
    // iterators are just fancy pointers, to be used together with STL containers, such as vector, list, etc.
    list<Element>::iterator it;

    /// @note usage of iterators
    //      (1) a.begin() is an iterator to the first element of a
    //      (2) a.end()   is an iterator BEYOND the last element of a
    //      (3) list iterators can be incremented, which means pointing to the next element
    //      (4) iterators are dereferenced like pointers: if `it` is an iterator, then `*it` is the element
    //      (5) you can access member variables/functions with arrow operators after iterators 
    // see the following usage
    for (it = a.begin(); it != a.end(); it++) {
        cout << (*it).m_val << " ";    // dereference and access member variable
        // cout << it->m_val << " ";   // same as above
    }
    cout << "\n";


    /// @note a good thing about iterators is that it can be used with container APIs
    it = a.begin();
    a.erase(it);    // erases the element that `it` points to (the first element)
    for (it = a.begin(); it != a.end(); it++) {
        cout << (*it).m_val << " ";
    }
    cout << "\n";

}