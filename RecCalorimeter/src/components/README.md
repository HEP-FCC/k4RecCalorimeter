| Feature              | FCCSW `CaloDIGI`                   | CMSSW `EcalDigiProducer`                   |
| -------------------- | ---------------------------------- | ------------------------------------------ |
| Framework            | Gaudi (FCCSW)                      | CMSSW (EDM)                                |
| Input                | `edm4hep::SimCalorimeterHit`       | `std::vector<PCaloHit>`                    |
| Output               | 1 energy + time + position per hit | 10 ADC samples per crystal (`EBDataFrame`) |
| Pulse shaping        | ❌ Not implemented                  | ✅ Realistic (uses `EcalShape`)             |
| Gain switching       | ❌ Fixed gain                       | ✅ Multi-gain ADC (MGPA)                    |
| ADC resolution       | ✅ Configurable (e.g., 12-bit)      | ✅ MGPA model with gain + saturation        |
| Time sampling        | ❌ Single smeared time              | ✅ 10 samples at 25 ns spacing              |
| Noise model          | ✅ Gaussian                         | ✅ Gaussian, channel-specific               |
| Pedestal             | ✅ Global                           | ✅ Channel-specific from DB                 |
| Geometry             | ✅ DD4hep (FCCSW)                   | ✅ `EcalGeometry` / `CaloGeometry`          |
| Dead channel support | Optional                           | ✅ From `EcalChannelStatus`                 |
