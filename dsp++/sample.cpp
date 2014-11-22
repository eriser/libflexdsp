#include <dsp++/snd/sample.h>
#include <dsp++/snd/convert.h>
#include <dsp++/snd/buffer.h>

#include <limits>
#include <cassert>
#include <cstring>

using namespace dsp::snd;

namespace {
template<class Int, class Res>
inline void read_pcm_as_float_(const dsp::snd::sample_layout& sl, const void* data, Res& res) {
	Int i;
	sl.read_pcm(data, i);
	res = sample_cast<Res>(i);
}

template<class Int, class Float>
inline void write_float_as_pcm_(const dsp::snd::sample_layout& sl, Float in, void* out) {
	Int i = sample_cast<Int>(in);
	sl.write_pcm(i, out);
}

template<class Res>
inline void read_pcm_signed_as_float(const dsp::snd::sample_layout& sl, const void* data, Res& res) {
	if (sl.container_bytes > sizeof(int32_t)) 
		read_pcm_as_float_<int64_t>(sl, data, res);
	else if (sl.container_bytes> sizeof(int16_t))
		read_pcm_as_float_<int32_t>(sl, data, res);
	else if (sl.container_bytes > sizeof(int8_t))
		read_pcm_as_float_<int16_t>(sl, data, res);
	else
		read_pcm_as_float_<int8_t>(sl, data, res);
}

template<class Float>
inline void write_float_as_pcm_signed(const dsp::snd::sample_layout& sl, Float f, void* data) {
	if (sl.container_bytes > sizeof(int32_t)) 
		write_float_as_pcm_<int64_t>(sl, f, data);
	else if (sl.container_bytes> sizeof(int16_t))
		write_float_as_pcm_<int32_t>(sl, f, data);
	else if (sl.container_bytes > sizeof(int8_t))
		write_float_as_pcm_<int16_t>(sl, f, data);
	else
		write_float_as_pcm_<int8_t>(sl, f, data);
}

template<class Res>
inline void read_pcm_unsigned_as_float(const dsp::snd::sample_layout& sl, const void* data, Res& res) {
	if (sl.container_bytes > sizeof(uint32_t)) 
		read_pcm_as_float_<uint64_t>(sl, data, res);
	else if (sl.container_bytes> sizeof(uint16_t))
		read_pcm_as_float_<uint32_t>(sl, data, res);
	else if (sl.container_bytes > sizeof(uint8_t))
		read_pcm_as_float_<uint16_t>(sl, data, res);
	else
		read_pcm_as_float_<uint8_t>(sl, data, res);
}

template<class Float>
inline void write_float_as_pcm_unsigned(const dsp::snd::sample_layout& sl, Float f, void* data) {
	if (sl.container_bytes > sizeof(uint32_t)) 
		write_float_as_pcm_<uint64_t>(sl, f, data);
	else if (sl.container_bytes> sizeof(uint16_t))
		write_float_as_pcm_<uint32_t>(sl, f, data);
	else if (sl.container_bytes > sizeof(uint8_t))
		write_float_as_pcm_<uint16_t>(sl, f, data);
	else
		write_float_as_pcm_<uint8_t>(sl, f, data);
}

template<class Float>
inline void read_ieee_as_float(const dsp::snd::sample_layout& sl, const void* data, Float& res) {
	if (4 == sl.container_bytes) {
		dsp::float32_t f;
		sl.read_ieee_float(data, f);
		res = static_cast<Float>(f);
	}
	else if (8 == sl.container_bytes) {
		dsp::float64_t f;
		sl.read_ieee_float(data, f);
		res = static_cast<Float>(f);
	}
	else 
		throw std::runtime_error("dsp::snd::sample_layout::read_float() IEEE 754 supports only 32 and 64-bit containers");
}

template<class Float>
inline void write_float_as_ieee(const dsp::snd::sample_layout& sl, Float in, void* data) {
	if (4 == sl.container_bytes) {
		dsp::float32_t f = sample_cast<dsp::float32_t>(in);
		sl.write_ieee_float(f, data);
	}
	else if (8 == sl.container_bytes) {
		dsp::float64_t f = sample_cast<dsp::float64_t>(in);
		sl.write_ieee_float(f, data);
	}
	else 
		throw std::runtime_error("dsp::snd::sample_layout::write_float() IEEE 754 supports only 32 and 64-bit containers");
}


template<class Float>
inline void read_sample_as_float(const dsp::snd::sample_layout& sl, const void* data, Float& res) {
	switch (sl.type) {
	case sample::type::ieee_float:
		read_ieee_as_float(sl, data, res);
		break;
	case sample::type::pcm_signed:
		read_pcm_signed_as_float(sl, data, res);
		break;
	case sample::type::pcm_unsigned:
		read_pcm_unsigned_as_float(sl, data, res);
		break;
	default:
		throw std::runtime_error("dsp::snd::sample_layout::read_float() unknown sample format");
	}
}

template<class Float>
inline void write_sample_as_float(const dsp::snd::sample_layout& sl, Float in, void* data) {
	switch (sl.type) {
	case sample::type::ieee_float:
		write_float_as_ieee(sl, in, data);
		break;
	case sample::type::pcm_signed:
		write_float_as_pcm_signed(sl, in, data);
		break;
	case sample::type::pcm_unsigned:
		write_float_as_pcm_unsigned(sl, in, data);
		break;
	default:
		throw std::runtime_error("dsp::snd::sample_layout::write_float() unknown sample format");
	}
}

}

#if !defined(DSP_ENDIAN_LITTLE) && !defined(DSP_ENDIAN_BIG)

namespace {
static dsp::snd::byte_order::label platform_test() {
	int16_t i = 1;
	int8_t buf[2];
	std::memcpy(buf, &i, 2);
	return ((buf[0] != 0) ? dsp::snd::byte_order::little_endian : dsp::snd::byte_order::big_endian);
}
}

const dsp::snd::byte_order::label dsp::snd::byte_order::platform = platform_test();

#endif

void dsp::snd::sample_layout::read_float(const void* data, float& out) const {
	read_sample_as_float(*this, data, out);
}

void dsp::snd::sample_layout::read_float(const void* data, double& out) const {
	read_sample_as_float(*this, data, out);
}

void dsp::snd::sample_layout::write_float(float in, void* data) const {
	write_sample_as_float(*this, in, data);
}

void dsp::snd::sample_layout::write_float(double in, void* data) const {
	write_sample_as_float(*this, in, data);
}

#if 0

static const char data[] = "\0\1\2\3\4\5\6\7";

static bool test() {
	dsp::snd::sample_layout s8_le(dsp::snd::sample::type::pcm_signed, 1, dsp::snd::byte_order::little_endian);
	long long res;
	s8_le.read_pcm_int_unnormalized(res, &data[0]);
	assert(res == 0);
	s8_le.read_pcm_int_unnormalized(res, &data[1]);
	assert(res == 1);

	dsp::snd::sample_layout s16_le(dsp::snd::sample::type::pcm_signed, 2, dsp::snd::byte_order::little_endian);
	s16_le.read_pcm_int_unnormalized(res, &data[0]);
	assert(res = 256);
	s16_le.read_pcm_int_unnormalized(res, &data[2]);
	assert(res = 3 * 256 + 2);

	dsp::snd::sample_layout s24_le(dsp::snd::sample::type::pcm_signed, 3, dsp::snd::byte_order::little_endian);
	s24_le.read_pcm_int_unnormalized(res, &data[0]);
	assert(res = 1 * 256 + 2 * 65536);
	s24_le.read_pcm_int_unnormalized(res, &data[3]);
	assert(res = 3 + 4 * 256 + 5 * 65536);

	return true;
}


static const bool t = test();

#endif 

void dsp::snd::convert_samples(const dsp::snd::sample_layout& sl_in, unsigned sample_stride_in, const void* in,
	const dsp::snd::sample_layout& sl_out, unsigned sample_stride_out, void* out, unsigned length) 
{
	const uint8_t* bi = static_cast<const uint8_t*>(in);
	uint8_t* bo = static_cast<uint8_t*>(out);
	//sample::type::label fmt;
	//if (sample::type::ieee_float == sl_in.type || sample::type::ieee_float == sl_out.type)
	//	fmt = sample::type::ieee_float;
	//else if (sample::type::pcm_signed == sl_in.type || sample::type::pcm_signed == sl_out.type)
	//	fmt = sample::type::pcm_signed;
	//else
	//	fmt = sample::type::pcm_unsigned;
	unsigned bytes = std::max(sl_in.container_bytes, sl_out.container_bytes);
	
	for (unsigned i = 0; i < length; ++i, bi += sample_stride_in, bo += sample_stride_out) 
	//switch (fmt) {
	//case sample::type::ieee_float:
		if (bytes > 4) {
			float64_t f;
			sl_in.read_float(bi, f);
			sl_out.write_float(f, bo);
		} 
		else {
			float32_t f;
			sl_in.read_float(bi, f);
			sl_out.write_float(f, bo);
		};
	//	break;
	//}
}

void dsp::snd::convert_samples(const dsp::snd::sample_layout& sl_in, const dsp::snd::buffer_layout& bl_in, const void* in,
	const dsp::snd::sample_layout& sl_out, const dsp::snd::buffer_layout& bl_out, void* out, unsigned length, unsigned channels)
{
	const uint8_t* bi = static_cast<const uint8_t*>(in);
	uint8_t* bo = static_cast<uint8_t*>(out);

	for (unsigned c = 0; c < channels; ++c, bi += bl_in.channel_stride, bo += bl_out.channel_stride) 
		convert_samples(sl_in, bl_in.sample_stride, bi, sl_out, bl_out.sample_stride, bo, length);
}
