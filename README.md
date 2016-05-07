# cvec

Cvec is a wrapper of SIMD intrinsics for handling 8bit-character array, providing load / store and alphabet conversions.

## Types and Macros

`cvec_t ` is a type which contains 32 elements of 8-bit-wide `char`. `cvec_t` can be loaded from and stored to memory with `_load_cvec(char const *ptr)` and `_store_cvec(char *ptr, cvec_t cv)` in unaligned way. Some utility (alphabet conversion and '\0' detection) macros are provided.


### Misc

* `_set_cvec`: broadcast an element to a vector
* `_zero_cvec`: clear a vector


### `\0` detection

* `_null_cvec`: returns non-zero if '\0' is found in the vector
* `_strlen_cvec`: strlen compatible macro. returns >= 32 if '\0' is not found in the vector.

### Conversion macros

* `_conv_a5_cvec`: ASCII -> 5bit conversion ({ 'A', 'a' } -> 0x01, ..., { 'Z', 'z' } -> 0x1b )
* `_conv_5a_cvec`: 5bit -> lowercase ASCII
* `_conv_5A_cvec`: 5bit -> uppercase ASCII
* `_conv_4t_cvec`: 4bit -> arbitrary 8-bit alphabet (with conversion table)
* `_conv_5t_cvec`: 5bit -> arbitrary 8-bit alphabet (with conversion table)


## License

MIT