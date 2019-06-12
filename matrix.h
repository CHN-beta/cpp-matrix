#pragma once
#include<cmath>
#include<iostream>
#include<string>
#include<random>
#include<memory.h>
#include<functional>
#include<ctime>

template<class num_t> class matrix
	// 实现简单的矩阵运算
{
private:
	num_t *val;
	// 存储矩阵的值，i行j列（从0记起） 的值为 value[i * size[1] + j]
	int *size = new int[2];
	// 矩阵的行数、列数
public:
	matrix();
	matrix(int, int);
	matrix(int*);
	// 传入行数、列数，矩阵值置为0
	matrix(int, int, num_t*);
	matrix(int, int, num_t);
	// 传入行数、列数，读入数据
	matrix(const matrix<num_t>&);
	~matrix();

	matrix<num_t> &resize(int, int);
	num_t &value(int, int) const;
	matrix<num_t> row(int) const;
	matrix<num_t> column(int) const;
	matrix<num_t> &set_random();
	// 将矩阵各元素值置为0到1的随机浮点数

	int n_row() const;
	int n_column() const;
	// 返回行列数
	int length() const;
	// 返回总元素数

	matrix<num_t>& operator=(const matrix<num_t>&);
	matrix<num_t>& operator=(num_t);
	bool operator==(const matrix<num_t>&) const;

	matrix<num_t>& operator+=(const matrix<num_t>&);
	matrix<num_t>& operator+=(num_t);
	matrix<num_t> operator+(const matrix<num_t>&) const;
	matrix<num_t> operator+(num_t) const;

	matrix<num_t>& operator*=(num_t);
	matrix<num_t> operator*(matrix<num_t>) const;
	matrix<num_t> operator*(num_t) const;

	matrix<num_t> operator-() const;
	matrix<num_t> operator-(const matrix<num_t>&) const;
	matrix<num_t>& operator-=(const matrix<num_t>&);
    matrix<num_t> operator-(const num_t);

	matrix<num_t>& operator/=(num_t);
	matrix<num_t> operator/(num_t) const;


	matrix<num_t> transposition() const;
	// 返回矩阵的转置
	num_t sum() const;
	// 返回矩阵各元素的和
	matrix<num_t> sum_row() const;
	matrix<num_t> ew_sign() const;
	// 返回矩阵各元素的符号

	matrix<num_t> ew_power_as_exponent(num_t) const;
	matrix<num_t> ew_exp() const;
	// 以数为底、以矩阵各元素为指数的乘方
	matrix<num_t> ew_power_as_base(num_t) const;
	// 以矩阵各元素为底、以数为指数的乘方
	matrix<num_t> ew_multiply(const matrix<num_t>&) const;
	matrix<num_t> ew_reciprocal() const;
	// 对矩阵各个元素求倒数
	matrix<num_t> ew_sin() const;
	matrix<num_t> ew_cos() const;
	matrix<num_t> ew_function(std::function<num_t(const num_t)>) const;
	// 对矩阵的各个元素使用自定义函数计算

	static matrix<num_t> linspace(num_t, num_t, int);
	// 传入上下限、个数，生成等差的列向量
	static matrix<num_t> random(int, int);
	// 生成指定大小的随机数矩阵
	static matrix<num_t> vander(matrix<num_t>, int);

	void print(std::ostream&) const;
};

template<class num_t> matrix<num_t>::matrix()
{
	size[0] = size[1] = 0;
	val = new num_t[0];
}
template<class num_t> matrix<num_t>::matrix(int n_row, int n_col)
{
	if (n_row < 0 || n_col < 0)
		throw("matirx size error");
	size[0] = n_row;
	size[1] = n_col;
	val = new num_t[length()];
	memset(val, 0, sizeof(num_t) * length());
}
template<class num_t> matrix<num_t>::matrix(int *size)
{
	if (size[0] < 0 || size[1] < 0)
		throw("matirx size error");
	memcpy(this->size, size, sizeof(int) * 2);
	val = new num_t[length()];
	memset(val, 0, sizeof(num_t) * length());
}
template<class num_t> matrix<num_t>::matrix(int n_row, int n_col, num_t *value)
{
	if (n_row < 0 || n_col < 0)
		throw("matirx size error");
	size[0] = n_row;
	size[1] = n_col;
	val = new num_t[length()];
	memcpy(val, value, sizeof(num_t) * length());
}
template<class num_t> matrix<num_t>::matrix(int n_row, int n_col, num_t value)
{
	if (n_row < 0 || n_col < 0)
		throw("matirx size error");
	size[0] = n_row;
	size[1] = n_col;
	val = new num_t[length()];
	for (int i = 0; i < length(); i++)
		val[i] = value;
}
template<class num_t> matrix<num_t>::matrix(const matrix<num_t> &ori)
{
	memcpy(size, ori.size, sizeof(int) * 2);
	val = new num_t[length()];
	memcpy(val, ori.val, sizeof(num_t) * length());
}
template<class num_t> matrix<num_t>::~matrix()
{
	delete[] val;
	delete[] size;
}

template<class num_t> matrix<num_t>& matrix<num_t>::resize(int i, int j)
{
	*this = matrix(i, j);
	return *this;
}
template<class num_t> num_t& matrix<num_t>::value(int i, int j) const
{
	if (i < 0 || i >= size[0] || j < 0 || j >= size[1])
		throw("matrix element index out of size");
	return val[i * size[1] + j];
}
template<class num_t> matrix<num_t> matrix<num_t>::row(int i) const
{
	if (i < 0 || i >= size[0])
		throw("matrix row index out of size");
	matrix<num_t> ans(1, size[1]);
	memcpy(ans.val, &value(i, 0), sizeof(num_t) * size[1]);
	return ans;
}
template<class num_t> matrix<num_t> matrix<num_t>::column(int i) const
{
	if (i < 0 || i >= size[1])
		throw("matrix column index out of size");
	matrix<num_t> ans(size[0], 1);
	for (int j = 0; j < size[0]; j++)
		ans.value(j, 0) = value(j, i);
	return ans;
}
template<class num_t> matrix<num_t>& matrix<num_t>::set_random()
{
	for (int i = 0; i < length(); i++)
		val[i] = ((num_t)rand()) / RAND_MAX * 2 - 1;
	return *this;
}

template<class num_t> int matrix<num_t>::n_row() const
{
	return size[0];
}
template<class num_t> int matrix<num_t>::n_column() const
{
	return size[1];
}
template<class num_t> int matrix<num_t>::length() const
{
	return size[0] * size[1];
}

template<class num_t> matrix<num_t>& matrix<num_t>::operator=(const matrix<num_t>& b)
{
	memcpy(size, b.size, sizeof(int) * 2);
	delete[] val;
	val = new num_t[length()];
	memcpy(val, b.val, sizeof(num_t) * length());
	return *this;
}
template<class num_t> matrix<num_t>& matrix<num_t>::operator=(num_t b)
{
	for (int i = 0; i < length(); i++)
		val[i] = b;
	return *this;
}
template<class num_t> bool matrix<num_t>::operator==(const matrix<num_t>& b) const
{
	if (memcmp(size, b.size, sizeof(int) * 2))
		return false;
	else
		return !memcmp(val, b.val, sizeof(num_t) * length());
}

template<class num_t> matrix<num_t>& matrix<num_t>::operator+=(const matrix<num_t>& b)
{
	if (memcmp(size, b.size, sizeof(int) * 2))
		throw("matrix size not match");
	for (int i = 0; i < length(); i++)
		val[i] += b.val[i];
	return *this;
}
template<class num_t> matrix<num_t>& matrix<num_t>::operator+=(num_t b)
{
	for (int i = 0; i < length(); i++)
		val[i] += b;
	return *this;
}
template<class num_t> matrix<num_t> matrix<num_t>::operator+(const matrix<num_t>& b) const
{
	matrix ans(*this);
	ans += b;
	return ans;
}
template<class num_t> matrix<num_t> matrix<num_t>::operator+(num_t b) const
{
	matrix ans(*this);
	ans += b;
	return ans;
}

template<class num_t> matrix<num_t>& matrix<num_t>::operator*=(num_t b)
{
	for (int i = 0; i < length(); i++)
		val[i] *= b;
	return *this;
}
template<class num_t> matrix<num_t> matrix<num_t>::operator*(matrix<num_t> b) const
{
	if (size[1] != b.size[0])
		throw("matrix size not match");
	matrix<num_t> ans(size[0], b.size[1]);
	for (int i = 0; i < size[0]; i++)
		for (int j = 0; j < b.size[1]; j++)
			for (int k = 0; k < size[1]; k++)
				ans.value(i, j) += value(i, k) * b.value(k, j);
	return ans;
}
template<class num_t> matrix<num_t> matrix<num_t>::operator*(num_t b) const
{
	matrix<num_t> ans(*this);
	ans *= b;
	return ans;
}

template<class num_t> matrix<num_t> matrix<num_t>::operator-() const
{
	return *this * (-1);
}
template<class num_t> matrix<num_t> matrix<num_t>::operator-(const matrix<num_t>& b) const
{
	return *this + (-b);
}
template<class num_t> matrix<num_t>& matrix<num_t>::operator-=(const matrix<num_t>& b)
{
	return *this += (-b);
}
template<class num_t> matrix<num_t> matrix<num_t>::operator-(const num_t b)
{
    return *this + (-b);
}

template<class num_t> matrix<num_t>& matrix<num_t>::operator/=(num_t b)
{
	for (int i = 0; i < length(); i++)
		val[i] /= b;
	return *this;
}
template<class num_t> matrix<num_t> matrix<num_t>::operator/(num_t b) const
{
	matrix<num_t> ans(*this);
	ans /= b;
	return ans;
}

template<class num_t> matrix<num_t> matrix<num_t>::transposition() const
{
	matrix<num_t> ans(size[1], size[0]);
	for (int i = 0; i < ans.size[0]; i++)
		for (int j = 0; j < ans.size[1]; j++)
			ans.value(i, j) = value(j, i);
	return ans;
}
template<class num_t> num_t matrix<num_t>::sum() const
{
	num_t ans = 0;
	for (int i = 0; i < length(); i++)
		ans += val[i];
	return ans;
}
template<class num_t> matrix<num_t> matrix<num_t>::sum_row() const
{
	matrix<num_t> ans(n_row(), 1, num_t(0));
	for (int i = 0; i < n_row(); i++)
		for (int j = 0; j < n_column(); j++)
			ans.value(i, 0) += value(i, j);
	return ans;
}
template<class num_t> matrix<num_t> matrix<num_t>::sign() const
{
    matrix<num_t> rtn(this->n_row(), this->n_column());
    for (int i = 0; i < this->length(); i++)
        if (this->val[i] > 0)
            rtn.val[i] = 1;
        else if (this->val[i] < 0)
            rtn.val[i] = -1;
        else
            rtn.val[i] = 0;
    return rtn;
}

template<class num_t> matrix<num_t> matrix<num_t>::ew_power_as_exponent(num_t a) const
{
	matrix<num_t> ans(size);
	for (int i = 0; i < ans.length(); i++)
		ans.val[i] = pow(a, val[i]);
	return ans;
}
template<class num_t> matrix<num_t> matrix<num_t>::ew_exp() const
{
	matrix<num_t> ans(size);
	for (int i = 0; i < length(); i++)
		ans.val[i] = exp(val[i]);
	return ans;
}
template<class num_t> matrix<num_t> matrix<num_t>::ew_power_as_base(num_t b) const
{
	matrix<num_t> ans(size);
	for (int i = 0; i < length(); i++)
		ans.val[i] = pow(val[i], b);
	return ans;
}
template<class num_t> matrix<num_t> matrix<num_t>::ew_multiply(const matrix<num_t>& b) const
{
	if (memcmp(size, b.size, sizeof(int) * 2))
		throw("matrix size not match");
	matrix<num_t> ans(size);
	for (int i = 0; i < length(); i++)
		ans.val[i] = val[i] * b.val[i];
	return ans;
}
template<class num_t> matrix<num_t> matrix<num_t>::ew_reciprocal() const
{
	matrix<num_t> ans(size);
	for (int i = 0; i < length(); i++)
		ans.val[i] = 1 / val[i];
	return ans;
}
template<class num_t> matrix<num_t> matrix<num_t>::ew_sin() const
{
	matrix<num_t> ans(size);
	for (int i = 0; i < length(); i++)
		ans.val[i] = sin(val[i]);
	return ans;
}
template<class num_t> matrix<num_t> matrix<num_t>::ew_cos() const
{
	matrix<num_t> ans(size);
	for (int i = 0; i < length(); i++)
		ans.val[i] = cos(val[i]);
	return ans;
}
template<class num_t> matrix<num_t> matrix<num_t>::ew_function(std::function<num_t(num_t)> f) const
{
	matrix<num_t> ans(size);
	for (int i = 0; i < length(); i++)
		ans.val[i] = f(val[i]);
	return ans;
}

template<class num_t> matrix<num_t> matrix<num_t>::linspace(num_t a, num_t b, int n)
{
	if (n < 0)
		throw("matrix element number error");
	matrix<num_t> ans(n, 1);
	for (int i = 0; i < n; i++)
		ans.val[i] = a + (b - a) / n * i;
	return ans;
}
template<class num_t> matrix<num_t> matrix<num_t>::random(int a, int b)
{
	matrix<num_t> rtn = matrix<num_t>(a, b);
	rtn.set_random();
	return rtn;
}
template<class num_t> matrix<num_t> matrix<num_t>::vander(matrix<num_t> a, int n)
{
	matrix<num_t> rtn = matrix<num_t>(a.length(), n);
	if (n > 0)
		for (int i = 0; i < a.length(); i++)
			rtn.value(i, 0) = 1;
	for (int i = 0; i < a.length(); i++)
		for (int j = 1; j < n; j++)
			rtn.value(i, j) = a.val[i] * rtn.value(i, j - 1);
	return rtn;
}

template<class num_t> void matrix<num_t>::print(std::ostream& out) const
{
	out << size[0] << ' ' << size[1] << '\n';
	for (int i = 0; i < size[0]; i++)
	{
		for (int j = 0; j < size[1]; j++)
			out << value(i, j) << " ";
		out << '\n';
	}
	out << std::endl;
}
