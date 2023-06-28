#pragma once

#include <array>
#include <iostream>
#include <iterator>
#include <type_traits>
#include <vector>

template <typename T>
class Deque {
  public:
    static const size_t BLOCK_SIZE = 32;

    struct Block {
        T* data;
        size_t begin;
        size_t end;

        Block();
        ~Block();

        T& operator[](size_t pos);
        const T& operator[](size_t pos) const;

        T& at(size_t pos);
        const T& at(size_t pos) const;

        void push_front(const T& val);
        void pop_front();

        void push_back(const T& val);
        void pop_back();
    };

    template <bool is_const>
    class common_iterator {
      public:
        using value_type = T;
        using reference = std::conditional_t<is_const, T const&, T&>;
        using pointer = std::conditional_t<is_const, const T*, T*>;
        using block_ptr = Block**;

        common_iterator() = default;

        common_iterator(const common_iterator& other) = default;

        common_iterator& operator=(const common_iterator& other) = default;

        reference operator*() const {
            return *ptr_;
        }

        pointer operator->() const {
            return ptr_;
        };

        reference operator[](ptrdiff_t x) const {
            ptrdiff_t diff = (ptr_ - (*block_)->data);
            ptrdiff_t shift;
            if (diff + x >= 0) {
                shift = (diff + x) / static_cast<ptrdiff_t>(BLOCK_SIZE);
            } else {
                shift = (diff + x - static_cast<ptrdiff_t>(BLOCK_SIZE) + 1) /
                        static_cast<ptrdiff_t>(BLOCK_SIZE);
            }
            return (block_ + shift)->data +
                   (diff + x + BLOCK_SIZE) % BLOCK_SIZE;
        }

        ptrdiff_t operator-(const common_iterator& other) const {
            return pos_ - other.pos_;
        }

        common_iterator& operator+=(ptrdiff_t x) {
            ptrdiff_t diff = (ptr_ - (*block_)->data);
            ptrdiff_t shift;
            if (diff + x >= 0) {
                shift = (diff + x) / static_cast<ptrdiff_t>(BLOCK_SIZE);
            } else {
                shift = (diff + x - static_cast<ptrdiff_t>(BLOCK_SIZE) + 1) /
                        static_cast<ptrdiff_t>(BLOCK_SIZE);
            }
            block_ += shift;
            ptr_ = (*block_)->data + (diff + x + BLOCK_SIZE) % BLOCK_SIZE;
            pos_ += x;
            return *this;
        }

        common_iterator& operator-=(ptrdiff_t x) {
            return operator+=(-1 * x);
        }
        common_iterator& operator++() {
            return operator+=(1);
        }
        common_iterator& operator--() {
            return operator-=(1);
        };

        common_iterator operator++(int) {
            auto copy = *this;
            operator+=(1);
            return copy;
        }

        common_iterator operator--(int) {
            auto copy = *this;
            operator-=(1);
            return copy;
        }

        bool operator==(const common_iterator& other) const {
            return block_ == other.block_ && ptr_ == other.ptr_;
        }

        bool operator!=(const common_iterator& other) const {
            return !(*this == other);
        }

        bool operator<(const common_iterator& other) const {
            return block_ < other.block_ ||
                   (block_ == other.block_ && ptr_ < other.ptr_);
        }

        bool operator>(const common_iterator& other) const {
            return other < *this;
        }

        bool operator<=(const common_iterator& other) const {
            return !(other < *this);
        }

        bool operator>=(const common_iterator& other) const {
            return !(other > *this);
        }

        friend common_iterator operator+(common_iterator lhs, ptrdiff_t rhs) {
            return lhs += rhs;
        }

        friend common_iterator operator+(ptrdiff_t lhs, common_iterator rhs) {
            return rhs += lhs;
        }

        common_iterator operator-(ptrdiff_t rhs) const {
            auto copy = *this;
            return copy -= rhs;
        }

        operator common_iterator<true>() const {
            return common_iterator<true>(ptr_, pos_, block_);
        }

      private:
        explicit common_iterator(pointer ptr, size_t pos, block_ptr block)
            : ptr_(ptr), pos_(pos), block_(block) {}

        pointer ptr_;
        size_t pos_;
        block_ptr block_;
        friend Deque;
    };

    using iterator = common_iterator<false>;
    using const_iterator = common_iterator<true>;
    using reverse_iterator = std::reverse_iterator<iterator>;
    using const_reverse_iterator = std::reverse_iterator<const_iterator>;

    Deque();
    Deque(size_t sz);
    Deque(size_t sz, const T& val);
    Deque(const Deque<T>& deq);
    ~Deque();

    Deque<T>& operator=(Deque<T> deq);

    T& operator[](size_t pos);
    const T& operator[](size_t pos) const;

    T& at(size_t pos);
    const T& at(size_t pos) const;

    iterator begin();
    const_iterator begin() const;
    std::reverse_iterator<iterator> rbegin();
    std::reverse_iterator<const_iterator> rbegin() const;

    const_iterator cbegin() const;
    std::reverse_iterator<const_iterator> crbegin() const;

    iterator end();
    const_iterator end() const;
    std::reverse_iterator<iterator> rend();
    std::reverse_iterator<const_iterator> rend() const;

    const_iterator cend() const;
    std::reverse_iterator<const_iterator> crend() const;

    void insert(iterator iter, const T& val);
    void erase(iterator iter);

    void push_back(const T& val);
    void pop_back();
    void push_front(const T& val);
    void pop_front();

    size_t size() const;
    void resize(size_t new_cap);

  private:
    void ClearBlocks();

    size_t blocks_begin_;
    size_t blocks_end_;
    size_t size_;
    size_t cap_;
    Block** blocks_;
};

template <typename T>
Deque<T>::Block::Block() : begin(0), end(0) {
    data = reinterpret_cast<T*>(new char[BLOCK_SIZE * sizeof(T)]);
}

template <typename T>
Deque<T>::Block::~Block() {
    for (size_t i = begin; i < end; ++i) {
        (data + i)->~T();
    }
    delete[] reinterpret_cast<char*>(data);
}

template <typename T>
T& Deque<T>::Block::operator[](size_t pos) {
    return data[pos];
}

template <typename T>
const T& Deque<T>::Block::operator[](size_t pos) const {
    return data[pos];
}

template <typename T>
T& Deque<T>::Block::at(size_t pos) {
    if (begin <= pos && pos < end) {
        return data[pos];
    }
    throw std::out_of_range("Error! Deque index out of range!");
}

template <typename T>
const T& Deque<T>::Block::at(size_t pos) const {
    if (begin <= pos && pos < end) {
        return data[pos];
    }
    throw std::out_of_range("Error! Deque index out of range!");
}

template <typename T>
void Deque<T>::Block::push_front(const T& val) {
    try {
        size_t new_begin = begin;
        size_t new_end = end;
        if (new_end == 0) {
            new_begin = BLOCK_SIZE;
            new_end = BLOCK_SIZE;
        }
        --new_begin;
        new (data + new_begin) T(val);
        begin = new_begin;
        end = new_end;
    } catch (...) {
        throw;
    }
}

template <typename T>
void Deque<T>::Block::pop_front() {
    (data + begin)->~T();
    ++begin;
}

template <typename T>
void Deque<T>::Block::push_back(const T& val) {
    size_t new_end = end;
    new (data + new_end) T(val);
    ++new_end;
    end = new_end;
}

template <typename T>
void Deque<T>::Block::pop_back() {
    (data + end - 1)->~T();
    if (end != 0) {
        --end;
    }
}

template <typename T>
Deque<T>::Deque() : blocks_begin_(1), blocks_end_(2), size_(0), cap_(2) {
    blocks_ = new Block*[2];
    blocks_[0] = new Block();
    blocks_[1] = new Block();
}

template <typename T>
Deque<T>::Deque(size_t sz)
    : blocks_begin_(0),
      blocks_end_((sz + BLOCK_SIZE) / BLOCK_SIZE),
      size_(sz),
      cap_((sz + BLOCK_SIZE) / BLOCK_SIZE + 1),
      blocks_(new Block*[cap_]) {
    for (size_t i = 0; i < cap_; ++i) {
        blocks_[i] = nullptr;
    }
    size_t i = blocks_begin_;
    try {
        for (; i <= blocks_end_; ++i) {
            blocks_[i] = new Block();
            blocks_[i]->end = (i != blocks_end_ - 1) ? BLOCK_SIZE : sz;
            sz -= BLOCK_SIZE;
            if constexpr (std::is_default_constructible_v<T>) {
                size_t k = 0;
                try {
                    for (; k < blocks_[i]->end; ++k) {
                        new (blocks_[i]->data + k) T();
                    }
                } catch (...) {
                    for (size_t j = 0; j < k; ++j) {
                        (blocks_[i]->data + j)->~T();
                    }
                    blocks_[i]->end = 0;
                    delete blocks_[i];
                    blocks_[i] = nullptr;
                    throw;
                }
            }
        }
    } catch (...) {
        ClearBlocks();
        throw;
    }
}

template <typename T>
Deque<T>::Deque(size_t sz, const T& val)
    : blocks_begin_(0),
      blocks_end_((sz + BLOCK_SIZE) / BLOCK_SIZE),
      size_(sz),
      cap_((sz + BLOCK_SIZE) / BLOCK_SIZE + 1),
      blocks_(new Block*[cap_]) {
    for (size_t i = 0; i < cap_; ++i) {
        blocks_[i] = nullptr;
    }

    size_t i = blocks_begin_;
    try {
        size_t cur = 0;
        for (; i <= blocks_end_; ++i) {
            blocks_[i] = new Block();
            size_t st = cur;
            for (; cur < std::min(size_, st + BLOCK_SIZE); ++cur) {
                blocks_[i]->push_back(val);
            }
        }
    } catch (...) {
        ClearBlocks();
        throw;
    }
}

template <typename T>
Deque<T>::Deque(const Deque<T>& deq)
    : blocks_begin_(deq.blocks_begin_),
      blocks_end_(deq.blocks_end_),
      size_(deq.size_),
      cap_(deq.cap_),
      blocks_(new Block*[deq.cap_]) {
    for (size_t i = 0; i < cap_; ++i) {
        blocks_[i] = nullptr;
    }

    size_t i = blocks_begin_;
    try {
        size_t cur = 0;
        for (; i < blocks_end_; ++i) {
            blocks_[i] = new Block();
            size_t st = cur;
            for (; cur < std::min(size_, st + BLOCK_SIZE); ++cur) {
                blocks_[i]->push_back(deq[cur]);
            }
        }
    } catch (...) {
        ClearBlocks();
        throw;
    }
}

template <typename T>
Deque<T>& Deque<T>::operator=(Deque<T> deq) {
    std::swap(blocks_begin_, deq.blocks_begin_);
    std::swap(blocks_end_, deq.blocks_end_);
    std::swap(size_, deq.size_);
    std::swap(cap_, deq.cap_);
    std::swap(blocks_, deq.blocks_);
    return *this;
}

template <typename T>
Deque<T>::~Deque() {
    ClearBlocks();
}

template <typename T>
T& Deque<T>::operator[](size_t pos) {
    pos += blocks_[blocks_begin_]->begin;
    size_t shift = pos / BLOCK_SIZE;
    return blocks_[blocks_begin_ + shift]->data[pos % BLOCK_SIZE];
}

template <typename T>
const T& Deque<T>::operator[](size_t pos) const {
    pos += blocks_[blocks_begin_]->begin;
    size_t shift = pos / BLOCK_SIZE;
    return blocks_[blocks_begin_ + shift]->data[pos % BLOCK_SIZE];
}

template <typename T>
T& Deque<T>::at(size_t pos) {
    if (size_ == 0) {
        throw std::out_of_range("Error! Deque index out of range!");
    }
    pos += blocks_[blocks_begin_]->begin;
    size_t shift = pos / BLOCK_SIZE;
    if (blocks_begin_ + shift >= blocks_end_) {
        throw std::out_of_range("Error! Deque index out of range!");
    }
    return blocks_[blocks_begin_ + shift]->at(pos % BLOCK_SIZE);
}

template <typename T>
const T& Deque<T>::at(size_t pos) const {
    if (size_ == 0) {
        throw std::out_of_range("Error! Deque index out of range!");
    }
    pos += blocks_[blocks_begin_]->begin;
    size_t shift = pos / BLOCK_SIZE;
    if (blocks_begin_ + shift >= blocks_end_) {
        throw std::out_of_range("Error! Deque index out of range!");
    }
    return blocks_[blocks_begin_ + shift]->at(pos % BLOCK_SIZE);
}

template <typename T>
typename Deque<T>::iterator Deque<T>::begin() {
    return iterator(
        &blocks_[blocks_begin_]->data[blocks_[blocks_begin_]->begin], 0,
        &blocks_[blocks_begin_]);
}

template <typename T>
typename Deque<T>::const_iterator Deque<T>::begin() const {
    return const_iterator(
        &blocks_[blocks_begin_]->data[blocks_[blocks_begin_]->begin], 0,
        &blocks_[blocks_begin_]);
}

template <typename T>
std::reverse_iterator<typename Deque<T>::iterator> Deque<T>::rbegin() {
    return std::reverse_iterator(end());
}

template <typename T>
std::reverse_iterator<typename Deque<T>::const_iterator> Deque<T>::rbegin()
    const {
    return std::reverse_iterator(end());
}

template <typename T>
typename Deque<T>::const_iterator Deque<T>::cbegin() const {
    return const_iterator(
        &blocks_[blocks_begin_]->data[blocks_[blocks_begin_]->begin], 0,
        &blocks_[blocks_begin_]);
}

template <typename T>
std::reverse_iterator<typename Deque<T>::const_iterator> Deque<T>::crbegin()
    const {
    return std::reverse_iterator(cend());
}

template <typename T>
typename Deque<T>::iterator Deque<T>::end() {
    return iterator(
        &blocks_[blocks_end_ - 1]->data[blocks_[blocks_end_ - 1]->end], size_,
        &blocks_[blocks_end_ - 1]);
}

template <typename T>
typename Deque<T>::const_iterator Deque<T>::end() const {
    return const_iterator(
        &blocks_[blocks_end_ - 1]->data[blocks_[blocks_end_ - 1]->end], size_,
        &blocks_[blocks_end_ - 1]);
}

template <typename T>
std::reverse_iterator<typename Deque<T>::iterator> Deque<T>::rend() {
    return std::reverse_iterator(begin());
}

template <typename T>
std::reverse_iterator<typename Deque<T>::const_iterator> Deque<T>::rend()
    const {
    return std::reverse_iterator(begin());
}

template <typename T>
typename Deque<T>::const_iterator Deque<T>::cend() const {
    return const_iterator(
        blocks_[blocks_end_ - 1]->data + blocks_[blocks_end_ - 1]->end, size_,
        &blocks_[blocks_end_ - 1]);
}

template <typename T>
std::reverse_iterator<typename Deque<T>::const_iterator> Deque<T>::crend()
    const {
    return std::reverse_iterator(cbegin());
}

template <typename T>
void Deque<T>::insert(iterator iter, const T& val) {
    size_t iter_pos = iter.pos_;
    push_back(val);
    iter = begin() + iter_pos;
    for (size_t i = size_ - 1; i > iter_pos; --i) {
        operator[](i) = operator[](i - 1);
    }
    operator[](iter_pos) = val;
}

template <typename T>
void Deque<T>::erase(iterator iter) {
    size_t iter_pos = iter.pos_;
    for (size_t i = iter_pos; i < size_ - 1; ++i) {
        operator[](i) = operator[](i + 1);
    }
    pop_back();
}

template <typename T>
void Deque<T>::push_back(const T& val) {
    blocks_[blocks_end_ - 1]->push_back(val);
    ++size_;
    if (blocks_[blocks_end_ - 1]->end == BLOCK_SIZE) {
        if (blocks_end_ == cap_) {
            resize(3 * cap_);
        } else {
            if (blocks_[blocks_end_] == nullptr) {
                blocks_[blocks_end_] = new Block();
            }
            ++blocks_end_;
        }
    }
}

template <typename T>
void Deque<T>::pop_back() {
    if (blocks_[blocks_end_ - 1]->end == 0) {
        delete blocks_[blocks_end_ - 1];
        blocks_[blocks_end_ - 1] = nullptr;
        --blocks_end_;
    }
    blocks_[blocks_end_ - 1]->pop_back();
    --size_;
}

template <typename T>
void Deque<T>::push_front(const T& val) {
    if (blocks_[blocks_begin_]->begin == 0) {
        if (blocks_begin_ == 0) {
            resize(3 * cap_);
        }
        --blocks_begin_;
        if (blocks_[blocks_begin_] == nullptr) {
            blocks_[blocks_begin_] = new Block();
        }
    }
    blocks_[blocks_begin_]->push_front(val);
    ++size_;
    if (blocks_[blocks_end_ - 1]->end == BLOCK_SIZE) {
        if (blocks_end_ == cap_) {
            resize(3 * cap_);
        } else {
            if (blocks_[blocks_end_] == nullptr) {
                blocks_[blocks_end_] = new Block();
            }
            ++blocks_end_;
        }
    }
}

template <typename T>
void Deque<T>::pop_front() {
    blocks_[blocks_begin_]->pop_front();
    --size_;
    if (blocks_[blocks_begin_]->begin == BLOCK_SIZE) {
        delete blocks_[blocks_begin_];
        blocks_[blocks_begin_] = nullptr;
        if (blocks_begin_ == blocks_end_ - 1) {
            blocks_[blocks_begin_] = new Block();
        } else {
            ++blocks_begin_;
        }
    }
}

template <typename T>
size_t Deque<T>::size() const {
    return size_;
}

template <typename T>
void Deque<T>::resize(size_t new_cap) {
    try {
        Block** new_blocks = new Block*[new_cap];
        size_t block_cnt = blocks_end_ - blocks_begin_;
        size_t new_blocks_begin = (new_cap - block_cnt) / 2;
        size_t new_blocks_end = (new_cap + block_cnt) / 2;
        if (new_blocks_end - new_blocks_begin != block_cnt) {
            ++new_blocks_end;
        }
        for (size_t i = 0; i < new_cap; ++i) {
            new_blocks[i] = nullptr;
        }
        size_t cur = new_blocks_begin;
        for (size_t i = blocks_begin_; i <= std::min(blocks_end_, cap_ - 1);
             ++i) {
            new_blocks[cur++] = std::move(blocks_[i]);
            blocks_[i] = nullptr;
        }
        for (size_t i = 0; i < cap_; ++i) {
            if (blocks_[i] != nullptr) {
                delete blocks_[i];
            }
        }
        delete[] blocks_;
        blocks_ = new_blocks;
        blocks_begin_ = new_blocks_begin;
        blocks_end_ = new_blocks_end;
        cap_ = new_cap;
        if (blocks_[blocks_end_ - 1]->end == BLOCK_SIZE) {
            blocks_[blocks_end_] = new Block();
            ++blocks_end_;
        }
    } catch (...) {
        throw;
    }
}

template <typename T>
void Deque<T>::ClearBlocks() {
    for (size_t i = 0; i < cap_; ++i) {
        if (blocks_[i] != nullptr) {
            delete blocks_[i];
        }
    }
    delete[] blocks_;
}
