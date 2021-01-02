#define _CRT_SECURE_NO_WARNINGS
#include<fstream>
#include<istream>
#define SR_NS_BEGIN(NS)     namespace NS {
#define SR_NS_END(NS)       }

#ifdef _XUTILITY_
#define ASSERT(expr, ...) if(!(expr)) __asm{ int 3 }
#else
#define ASSERT(expr, ...) 
#endif // DEBUG

SR_NS_BEGIN(euns)
SR_NS_END(euns)

//current if multiplication calculation exceeds defalutSize digits, program dead because memory allocation.

#include<memory>
#include<string>
#include<istream>
#include<ostream>
#include<iostream>
#include<string.h>
#include<complex>
#include<vector>
#include<assert.h>
#include<algorithm>

typedef std::complex<long double> base;
constexpr double PI = 3.14159265358979323846;


class bigint {
public:
	bigint();
	bigint(const char* const a);
	bigint(const std::string& a);

	bigint(const std::int8_t a); // 8 bit char
	bigint(const std::int16_t a); // 16bit short
	bigint(const std::int32_t a); //32bit int
	bigint(const std::int64_t a); // 64bit long long int
	bigint(const std::uint64_t a); // 64bit std::uint64_t

	//copy 
	bigint& operator=(const bigint& a);
	bigint& operator=(const std::string& a);
	bigint& operator=(const std::int64_t a);
	bigint& operator=(const std::uint64_t a);

	bigint(const bigint& a);
	//move
	/*bigint& operator=(bigint&& a) noexcept;
	bigint(bigint&& a) noexcept;*/
	~bigint();

	const bigint MultiplyFFT(const bigint& a) const;
	const bigint operator+(const bigint& a) const;
	const bigint operator-(const bigint& a) const;
	const bigint operator-() const; // -this 반환
	const bigint operator*(const bigint& a) const;// using fft
	const bigint operator/(const bigint& a) const;
	const bigint operator%(const bigint& a) const;
	bigint& operator+=(const bigint& a);
	bigint& operator-=(const bigint& a);
	bigint& operator*=(const bigint& a);
	bigint& operator/=(const bigint& a);
	bigint& operator%=(const bigint& a);

	/*operator std::int8_t() const;
	operator std::int64_t() const;
	operator std::int32_t() const;
	operator std::int16_t() const;*/

	bigint& operator++();
	bigint operator++(int);
	bigint& operator--();
	bigint operator--(int);

	bool operator==(const std::int64_t a) const;
	bool operator!=(const std::uint64_t  a) const;
	bool operator==(const bigint& a) const;
	bool operator!=(const bigint& a) const;
	bool operator>(const bigint& a) const;
	bool operator<(const bigint& a) const;
	bool operator>=(const bigint& a) const;
	bool operator<=(const bigint& a) const;

	std::string ToString();
	friend std::ostream& operator<<(std::ostream& os, const bigint& a);
	friend std::istream& operator>>(std::istream& os, bigint& a);
private:
	void fft(std::vector<base>& inOutF, bool invert) const;
	const bigint addition(const bigint& a) const;
	const bigint subtraction(const bigint& a) const;
private:
	static constexpr int mDEFULT_SIZE = 1020;
	static constexpr int mDIVIDE_NUM = 1000'000'000;
	static constexpr int mCARRYGE_NUM = 9;

	std::vector<std::int32_t> mString; // default 
	int mCapacity; //array size
	int mSize; //array max index
	bool mSign;//+ : true - : false 
};


bigint::bigint()
	: bigint(0)
{}

bigint::bigint(const std::string& a)
	: bigint(a.c_str())
{}


bigint::bigint(const char* a)
	: mCapacity(mDEFULT_SIZE),
	mSize(0),
	mSign(true)
{
	int numLength = 0;
	int havingSign = 0;
	if (*a == '-')
	{
		mSign = false;
		++havingSign;
	}
	else if (*a == '+')
	{
		++havingSign;
	}

	char reg = 0;
	while ((reg = *(a + havingSign + numLength)) != '\0')
	{
		if (('0' > reg || '9' < reg ))//&& reg != '\''
		{
			
			std::cout << "not number" << std::endl;
			ASSERT(false);
			return;
		}
		++numLength;
	}

	while (mCapacity * mCARRYGE_NUM <= numLength)
	{
		mCapacity *= 2;
	}
	
	mString.resize(mCapacity);
	mString[0] = 0;

	for (int j = 0; j < numLength; j += mCARRYGE_NUM)
	{
		int carry = numLength - j - mCARRYGE_NUM > 0 ? numLength - j - mCARRYGE_NUM : 0;
		int  y = numLength - j;
		for (int i = 0; i < mCARRYGE_NUM && y != 0; i++, y--)
		{
			int num = *(a + havingSign + carry + i) - '0';


			ASSERT(num >= 0 && num <= 9);
			mString[mSize] = mString[mSize] * 10 + num;
		}
		mSize++;
	}

}


bigint::bigint(const std::int8_t a) : bigint((std::int64_t)a)
{}

bigint::bigint(const std::int16_t a) : bigint((std::int64_t)a)
{}

bigint::bigint(const std::int32_t a) : bigint((std::int64_t)a)
{}

bigint::bigint(const std::int64_t a) // 64bit long long int
	: mSign(a >= 0),
	mSize(0),
	mCapacity(mDEFULT_SIZE)
{
	mString.resize(mCapacity);
	
	auto b = a >= 0 ? a : -a;
	ASSERT(b >= 0);
	if (a == 0)
	{
		mString[mSize] = 0;
		mSize++;
	}

	while (b)
	{
		mString[mSize] = b % mDIVIDE_NUM;
		b /= mDIVIDE_NUM;
		mSize++;
	}
}

bigint::bigint(const std::uint64_t a)
	: mSign(true),
	mSize(0),
	mCapacity(mDEFULT_SIZE)
{
	mString.resize(mDEFULT_SIZE);
	auto b = a;
	while (b)
	{
		mString[mSize] = b % mDIVIDE_NUM;
		b /= mDIVIDE_NUM;
		mSize++;
	}

}

// copy constructor
bigint::bigint(const bigint& a)
	: mCapacity(a.mCapacity),
	mSign(a.mSign),
	mSize(a.mSize),
	mString(a.mString)
{
	//ASSERT(a != NULL);
	//ASSERT(&a != NULL);
		//std::cout << "using copy constructor" << std::endl;
	//mString = std::make_unique<std::int32_t[]>(mCapacity);
	//memcpy(mString.get(), a.mString.get(), mCapacity*sizeof(std::int32_t));
	/*for (int i = 0; i < a.mCapacity; i++)
	{
		mString[i] = a.mString[i];
	}*/
}

bigint::~bigint()
{
}


bigint& bigint::operator=(const bigint& a)
{
	//return bigint(a);
	//std::cout << "nomal assignment operator" << '\n';
	if (&a == this)
	{
		return *this;
	}

	//if (mCapacity != a.mCapacity)
	//{
	//	mCapacity = a.mCapacity;
	//	mString.reset();
	//	mString = std::make_unique<std::int32_t[]>(mCapacity);
	//	for (int i = 0; i < mCapacity; i++)
	//	{
	//		mString[i] = 0;
	//	}
	//}
	//memcpy(mString.get(), a.mString.get(), mCapacity);
	mCapacity = a.mCapacity;
	mSign = a.mSign;
	mSize = a.mSize;
	mString = a.mString;
	return *this;
}

bigint& bigint::operator=(std::int64_t a)
{
	return (*this = bigint(a));
}

bigint& bigint::operator=(std::uint64_t a)
{
	return (*this = bigint(a));
}

bigint& bigint::operator=(const std::string& a)
{
	return (*this = bigint(a));
}

//
//bigint& bigint::operator=(bigint&& a) noexcept
//{
//	//	std::cout << "move assignment operator" << std::endl;
//	mString = std::move(a.mString);
//	mSign = (a.mSign);
//	mSize = (a.mSize);
//	mCapacity = (a.mCapacity);
//	a.mString = nullptr;
//	return *this;
//}
//
//bigint::bigint(bigint&& a)noexcept 
//	: mSign(a.mSign), 
//	mSize(a.mSize),
//	mCapacity(a.mCapacity)
//{
//	//	std::cout << "move constructor" << std::endl;
//	mString = std::move(a.mString);
//}

// only this procedure sign same
const bigint bigint::operator+(const bigint& a) const
{
	if (a.mSign != mSign)
	{
		return this->subtraction(a);
	}
	return this->addition(a);
}

// only this procedure sign deffernt
const bigint bigint::operator-(const bigint& a) const
{
	if (a.mSign != mSign)
	{
		return this->addition(a);
	}
	return this->subtraction(a);
}

const bigint bigint::addition(const bigint& a) const
{
	const bigint& tmp = *this > a ? (*this) : a;

	bigint returnValue(tmp);
	returnValue.mSign = mSign; // must this.mSign cuase operator-
	if (returnValue.mCapacity == returnValue.mSize)
	{
		returnValue.mCapacity *= 2;
		returnValue.mString.resize(returnValue.mCapacity);
	}

	int mmsize = a.mSize > mSize ? mSize : a.mSize;
	for (int i = 0; i < mmsize; i++)
	{
		returnValue.mString[i] = mString[i] + a.mString[i];
	}

	//ceiling
	for (int i = 0; i < returnValue.mSize; i++)
	{
		returnValue.mString[i + 1] += returnValue.mString[i] / mDIVIDE_NUM;
		returnValue.mString[i] = returnValue.mString[i] % mDIVIDE_NUM;
	}
	//ASSERT(returnValue.mSize != returnValue.mCapacity);
	if (returnValue.mString[returnValue.mSize] != 0)
	{
		returnValue.mSize++;
	}

	ASSERT(returnValue.mString[returnValue.mSize - 1] < mDIVIDE_NUM);
	return returnValue;
}


const bigint bigint::subtraction(const bigint& a) const
{
	const bigint* bbig = nullptr;
	const bigint* ssmall = nullptr;
	bool isReturnSign = true;

	//find abs bigger
	if (this->mSize == a.mSize)
	{
		for (int i = mSize - 1; i >= 0; i--)
		{
			if (mString[i] != a.mString[i])
			{
				if (mString[i] > a.mString[i])
				{
					bbig = this;
					ssmall = &a;
					break;
				}
				if (mString[i] < a.mString[i])
				{
					bbig = &a;
					ssmall = this;
					break;
				}
			}
		}
		if (bbig == nullptr && ssmall == nullptr)
		{
			// bbig == ssmall
			bbig = this;
			ssmall = &a;
		}

	}
	else
	{
		bbig = this->mSize > a.mSize ? this : &a;
		ssmall = this->mSize <= a.mSize ? this : &a;
	}

	ASSERT(bbig != nullptr && ssmall != nullptr);
	const bigint& big = *bbig;
	const bigint& small = *ssmall;

	bigint returnValue(big);

	if (big.mSign == small.mSign && bbig == &a)
	{
		returnValue.mSign = !a.mSign;
	}

	int carry = 0;
	for (int i = 0; i < small.mSize; i++)
	{
		returnValue.mString[i] = returnValue.mString[i] - small.mString[i] - carry;
		if (returnValue.mString[i] < 0)
		{
			returnValue.mString[i] += mDIVIDE_NUM;
			carry = 1;
		}
		else
		{
			carry = 0;
		}
	}

	for (int i = small.mSize; i < returnValue.mSize && carry != 0; i++)
	{
		returnValue.mString[i] = returnValue.mString[i] - carry;
		if (returnValue.mString[i] < 0)
		{
			returnValue.mString[i] += mDIVIDE_NUM;
			carry = 1;
		}
		else
		{
			carry = 0;
		}
	}

	ASSERT(carry == 0);

	int notZeroIndex = 0;
	for (int i = 0; i < returnValue.mSize; i++)
	{
		if (returnValue.mString[i] != 0)
		{
			notZeroIndex = i;
		}
	}
	returnValue.mSize = notZeroIndex + 1;
	if (returnValue.mSize == 1 && returnValue.mString[0] == 0)
	{
		returnValue.mSign = true;
	}

	return returnValue;
}

const bigint bigint::operator-() const
{
	//	std::cout << "test ";
	bigint returnValue = *this;
	if (returnValue != 0)
	{
		returnValue.mSign = !mSign;
	}
	return returnValue;
}

//TODO memroy 
const bigint bigint::operator*(const bigint& a) const
{
	if (a.mSize == 0 || mSize == 0)
	{
		return bigint(0);
	}

	bigint returnValue;
	returnValue.mCapacity = mCapacity > a.mCapacity ? mCapacity * 2 : a.mCapacity * 2;
	returnValue.mString.resize(returnValue.mCapacity);
	returnValue.mSign = (mSign == a.mSign);

	for (int i = 0; i < mSize; i++)
	{
		for (int j = 0; j < a.mSize; j++)
		{
			std::int64_t tmp = (std::int64_t)mString[i] * a.mString[j];
			returnValue.mString[i + j] = tmp % mDIVIDE_NUM;
			returnValue.mString[i + j + 1] += static_cast<std::int32_t>(tmp / mDIVIDE_NUM);
		}
	}
	int notZeroIndex = 0;
	for (int i = 0; i < returnValue.mCapacity; i++)
	{
		if (returnValue.mString[i] != 0)
		{
			notZeroIndex = i;
		}
	}
	returnValue.mSize = notZeroIndex + 1;
	if (returnValue.mSize == 1 && returnValue.mString[0] == 0)
	{
		returnValue.mSign = true;
	}

	return returnValue;
}

//Using fft Schönhage–Strassen algorithm Cooley-Tukey algorithm
// cause memory exception
const bigint bigint::MultiplyFFT(const bigint& a) const
{
	if (a.mSize == 0 || mSize == 0)
	{
		return bigint(0);
	}

	std::vector <base> fa(mString.begin(), mString.end() + mSize);
	std::vector <base>	fb(a.mString.begin(), a.mString.end() + a.mSize);
	int n = 1;
	while (n < mSize + a.mSize)
	{
		n <<= 1;
	}
	fa.resize(n);
	fb.resize(n);
	fft(fa, false);
	fft(fb, false);
	for (int i = 0; i < n; i++)
	{
		fa[i] *= fb[i];
	}
	fft(fa, true);

	std::vector<int> res;
	res.resize(n);
	for (int i = 0; i < n; i++)
	{
		res[i] = (int)round(fa[i].real());
	}

	while (!res.empty() && res.back() == 0)
	{
		res.pop_back();
	}

	for (int i = 0; i < (int)res.size() - 1; i++)
	{
		if (i == res.size() - 1)
		{
			if (res[i] >= 10)
			{

				res.push_back(res[res.size() - 1] / 10);
				res[res.size() - 2] %= 10;
			}
		}
		else {
			res[i + 1] += res[i] / 10;
			res[i] %= 10;
		}
	}
	bigint returnValue;
	returnValue.mSign = (a.mSign == mSign);
	if (returnValue.mCapacity < n)
	{
		returnValue.mString.resize(n);
		
	}

	for (int i = 0; i < (int)res.size(); i++)
	{
		returnValue.mString[i] = (char)res[i];
	}

	returnValue.mSize = res.size();

	return returnValue;
}

void bigint::fft(std::vector<base>& inOutF, bool invert) const
{
	const int n = inOutF.size();
	for (int i = 1, j = 0; i < n; i++)
	{
		int bit = n >> 1;
		for (; j >= bit; bit >>= 1)
		{
			j -= bit;
		}
		j += bit;
		if (i < j)
		{
			swap(inOutF[i], inOutF[j]);
		}
	}
	for (int len = 2; len <= n; len <<= 1)
	{
		double ang = 2 * PI / len * (invert ? -1 : 1);
		base wlen(cos(ang), sin(ang));
		for (int i = 0; i < n; i += len)
		{
			base w(1);
			for (int j = 0; j < len / 2; j++)
			{
				base u = inOutF[i + j];
				base v = inOutF[i + j + len / 2] * w;
				inOutF[i + j] = u + v;
				inOutF[i + j + len / 2] = u - v;
				w *= wlen;
			}
		}
	}

	if (invert)
	{
		for (int i = 0; i < n; i++)
		{
			inOutF[i] /= n;
		}
	}
}

//TODO asdfasdf
const bigint bigint::operator/(const bigint& a) const
{
	if (*this < a)
	{
		return bigint(0);
	}
	if (*this == a)
	{
		return bigint(1);
	}
	ASSERT(a < *this);

	bigint digit;
	int index = mSize - 1;
	int divideNum = mDIVIDE_NUM / 10;
	int sum = mString[index];

	while (digit < a)
	{
		if (divideNum == 0)
		{
			// ddigit == 0
			divideNum = mDIVIDE_NUM / 10;
			index--;
			sum = mString[index];
		}

		digit = digit * 10 + sum / divideNum;;

		sum %= divideNum;
		divideNum /= 10;

	}

	ASSERT(digit >= a);
	ASSERT(index >= 0);
	bigint returnValue;
	while (index >= 0)
	{

		for (int j = 9; j > 0; j--)
		{
			bigint reg = a * j;
			if (digit > reg)
			{
				digit -= reg;
				returnValue += j;
				break;
			}
		}

		if (divideNum == 0)
		{
			divideNum = mDIVIDE_NUM / 10;
			index--;
			if (index < 0)
			{
				break;
			}
			sum = mString[index];
		}

		digit = digit * 10 + sum / divideNum;
		sum %= divideNum;
		divideNum /= 10;
		returnValue *= 10;
	}
	return returnValue;
}

const bigint bigint::operator%(const bigint& a) const
{
	return (*this - (*this / a) * a);
}

bigint& bigint::operator+=(const bigint& a)
{
	return (*this = *this + a);
}

bigint& bigint::operator-=(const bigint& a)
{
	return (*this = *this - a);
}

bigint& bigint::operator*=(const bigint& a)
{
	return (*this = *this * a);
}


bigint& bigint::operator/=(const bigint& a)
{
	return (*this = *this / a);
}


bigint& bigint::operator%=(const bigint& a)
{
	return (*this = *this % a);
}

bigint& bigint::operator++()
{
	return (*this += 1);
}


bigint bigint::operator++(int)
{
	bigint returnValue = *this;
	++(*this);
	return returnValue;
}

bigint& bigint::operator--()
{
	return (*this -= 1);
}
bigint bigint::operator--(int)
{
	bigint returnValue = *this;
	--(*this);
	return returnValue;

}

bool bigint::operator==(std::int64_t a) const
{
	return *this == bigint(a);
}

bool bigint::operator!=(std::uint64_t  a) const
{
	return *this != bigint(a);
}


bool bigint::operator==(const bigint& a) const
{
	if (mSign != a.mSign)
	{
		return false;
	}
	if (mSize != a.mSize)
	{
		return false;
	}

	for (int i = mSize - 1; i >= 0; i--)
	{
		if (mString[i] != a.mString[i])
		{
			return false;
		}
	}
	return true;

}
bool bigint::operator!=(const bigint& a) const
{
	return !(*this == a);
}


bool bigint::operator>(const bigint& a) const
{
	if (mSign ^ a.mSign)
	{
		return mSign;
	}

	if (mSize != a.mSize)
	{
		return mSize > a.mSize;
	}

	for (int i = mSize - 1; i >= 0; i--)
	{
		if (mString[i] != a.mString[i])
		{
			return mString[i] > a.mString[i];
		}
	}
	return false;
}
bool bigint::operator<(const bigint& a) const
{
	if (mSign ^ a.mSign)
	{
		return !mSign;
	}

	if (mSize != a.mSize)
	{
		return mSize < a.mSize;
	}

	for (int i = mSize - 1; i >= 0; i--)
	{
		if (mString[i] != a.mString[i])
		{
			return mString[i] < a.mString[i];
		}
	}
	return false;
}

bool bigint::operator>=(const bigint& a) const
{
	return !(*this < a);
}
bool bigint::operator<=(const bigint& a) const
{
	return !(*this > a);
}
std::string bigint::ToString()
{
	std::string returnValue;
	
	if (mSign == false)
	{
		returnValue += '-';
	}
	returnValue += std::to_string(mString[mSize - 1]);
	for (int i = 0; i < mSize - 1; i++)
	{
		if (mString[mSize - 2 - i] < mDIVIDE_NUM / 10)
		{
			// (int i = 0; i < a.mCARRYGE_NUM; i++)
			int cnt = mDIVIDE_NUM / 10;
			while (mString[mSize - 2 - i] < cnt)
			{
				returnValue += '0';
				cnt /= 10;
			}//max 8 while
			if (mString[mSize - 2 - i] != 0)
			{
			returnValue += std::to_string(mString[mSize - 2 - i]);
			}
		}
		else
		{
			returnValue += std::to_string(mString[mSize - 2 - i]);
		}
	}

	return returnValue;
}

std::ostream& operator<<(std::ostream& os, const bigint& a)
{

	const int len = a.mSize;
	ASSERT((len != 1 && a.mString[len - 1] != 0) || len == 1);

	if (a.mSign == false)
	{
		os << '-';
	}
	int k = a.mString[len - 1];
	os << k;
	for (int i = 0; i < len - 1; i++)
	{
		if (a.mString[len - 2 - i] < a.mDIVIDE_NUM/10)
		{
			// (int i = 0; i < a.mCARRYGE_NUM; i++)
			int cnt = a.mDIVIDE_NUM / 10;
			while(a.mString[len - 2 - i] < cnt)
			{
				os << 0;
				cnt /= 10;
			}//max 8 while
			if (a.mString[len - 2 - i] != 0)
			{
				os << a.mString[len - 2 - i];
			}
			
		}
		else
		{
			os << a.mString[len - 2 - i];
		}
	}
	return os;
}
std::istream& operator>>(std::istream& os, bigint& a)
{
	std::string tmp;
	os >> tmp;
	a = tmp;
	return os;
}

char arr[100003];

int main()
{
	using namespace std;
#ifdef _XUTILITY_
	freopen("input.txt", "r", stdin);
//	freopen("output.txt", "w", stdout);
#endif // DEBUG
	/*for (int j = -100; j < 100; j++)
	{
		for (int i = -100; i < 100; i++)
		{*/

	/*std::ofstream out1("output1.txt");
	std::ofstream out2("output2.txt");
	*/
	/*cin >> arr;
	bigint a(arr);

	cout << a % 20000303;*/


	cin >> arr;
	bigint a(arr);
	cin >> arr;
	bigint b(arr);
	cout << a + b << '\n';
	cout << a - b << '\n';
	cout << a * b << '\n';

		/*cout << "-616512316543123164849124248688530900100237520000000000" << '\n';
		cout << c.ToString().c_str();
		if ( !strcmp(c.ToString().c_str() , "-616512316543123164849124248688530900100237520000000000"))
		{
			cout << 3;
		}*/



		//out2 << d + e << ' ';
	//}
	//out1.close();
	//out2.close();

	//ifstream in1("output1.txt");
	//ifstream in2("output2.txt");
	//d = -10;
	//while ((d = d + 1000000) < -e)
	//{
	//	long long t,w;
	//	in1 >> t;
	//	in2 >> w;
	//	if (t != w)
	//	{
	//		cout << w << '\n';
	//	}
	//}
	
		// 616512316543123166082148881774777231031070610811281290
	return 0;
}
