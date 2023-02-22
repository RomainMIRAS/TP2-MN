#include "complexe.h"

complexe_float_t add_complexe_float (const complexe_float_t c1, const complexe_float_t c2)
{
  complexe_float_t r ;

  r.real = c1.real + c2.real ;
  r.imaginary = c1.imaginary + c2.imaginary ;
  
  return r ;
}

complexe_double_t add_complexe_double (const complexe_double_t c1, const complexe_double_t c2)
{
  complexe_double_t r ;

  r.real = c1.real + c2.real ;
  r.imaginary = c1.imaginary + c2.imaginary ;
  
  return r ;
}

//6 opérations
complexe_float_t mult_complexe_float (const complexe_float_t c1, const complexe_float_t c2)
{
  complexe_float_t r ;

  r.real = c1.real * c2.real - c1.imaginary * c2.imaginary;

  r.imaginary = c1.imaginary * c2.real + c1.real * c2.imaginary;
  
  return r ;
}

complexe_double_t mult_complexe_double (const complexe_double_t c1, const complexe_double_t c2)
  {
  complexe_double_t r ;

  r.real = c1.real * c2.real - c1.imaginary * c2.imaginary;

  r.imaginary = c1.imaginary * c2.real + c1.real * c2.imaginary;
  
  return r ;
}
  

complexe_float_t div_complexe_float (const complexe_float_t c1, const complexe_float_t c2)
{
  complexe_float_t r ;

  complexe_float_t inverseBot;
  inverseBot.imaginary = -c2.imaginary;
  inverseBot.real = c2.real;

  complexe_float_t top = mult_complexe_float(c1, inverseBot);
  complexe_float_t bot = mult_complexe_float(c2, inverseBot);

  r.real = top.real/bot.real;
  r.imaginary = bot.imaginary/bot.real;

  
  return r ;
}

complexe_double_t div_complexe_double (const complexe_double_t c1, const complexe_double_t c2)
{
  complexe_double_t r ;

  complexe_double_t inverseBot;
  inverseBot.imaginary = -c2.imaginary;
  inverseBot.real = c2.real;

  complexe_double_t top = mult_complexe_double(c1, inverseBot);
  complexe_double_t bot = mult_complexe_double(c2, inverseBot);

  r.real = top.real/bot.real;
  r.imaginary = bot.imaginary/bot.real;
  
  return r ;
}
