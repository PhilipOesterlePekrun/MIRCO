namespace MIRCO
{
  // Shape factors (See section 3.3 of https://doi.org/10.1007/s00466-019-01791-3)
  // These are the shape factors to calculate the elastic compliance correction of the micro-scale
  // contact constitutive law for various N (number of elements/cells on the side length). These
  // were calculated with the flatMirco utility.


  constexpr double shape_factors_pressure[] = {-1, 1.1221997046783601, 1.0068605251532485,
      0.9613892379176019, 0.9385553774805636, 0.9247153424324350, 0.9154315222403194,
      0.9087685392496161, 0.9037519877715218, 0.8998375318806964, 0.8966971707662920,
      0.8941214388008705, 0.8919702710377435, 0.8901464435834811, 0.8885803476015518,
      0.8872208283196537, 0.8860294454171725, 0.8849767510419421, 0.8840398035738947,
      0.8832004634308468, 0.8824441981346327, 0.8817592274578133, 0.8811359008070114,
      0.8805662363985260, 0.8800435751941625, 0.8795623175713739, 0.8791177205270474,
      0.8787057397731856};



  // Originally present shape factors for integer resolutions of 1 to 8:
  /*
  // The following pressure based constants are calculated by solving a flat indentor problem using
  // the pressure based Green function described in Pohrt and Li (2014).
  // http://dx.doi.org/10.1134/s1029959914040109
  const std::map<int, double> shape_factors_pressure{{1, 0.961389237917602}, {2, 0.924715342432435},
      {3, 0.899837531880697}, {4, 0.884976751041942}, {5, 0.876753783192863},
      {6, 0.872397956576882}, {7, 0.8701463093314326}, {8, 0.8689982669426167}};

  // The following force based constants are taken from Table 1 of Bonari et al. (2020).
  // https://doi.org/10.1007/s00466-019-01791-3
  const std::map<int, double> shape_factors_force{{1, 0.778958541513360}, {2, 0.805513388666376},
      {3, 0.826126871395416}, {4, 0.841369158110513}, {5, 0.851733020725652},
      {6, 0.858342234203154}, {7, 0.862368243479785}, {8, 0.864741597831785}};
  */
}  // namespace MIRCO
