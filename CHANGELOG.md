### Changelog

All notable changes to this project will be documented in this file. Dates are displayed in UTC.

#### [v1.4.3](https://github.com/emit-sds/emit-sds-l1b/compare/v1.4.2...v1.4.3)

> 2 April 2025

* New CPM calibration scripts and results by @davidraythompson in https://github.com/emit-sds/emit-sds-l1b/pull/21
* removed CPM by @davidraythompson in https://github.com/emit-sds/emit-sds-l1b/pull/23
* New OSF interpolation strategy by @davidraythompson in https://github.com/emit-sds/emit-sds-l1b/pull/24
* Calibration update  by @davidraythompson in https://github.com/emit-sds/emit-sds-l1b/pull/26
* Merge main into develop by @winstonolson in https://github.com/emit-sds/emit-sds-l1b/pull/27

#### [v1.4.2](https://github.com/emit-sds/emit-sds-l1b/compare/v1.4.1...v1.4.2)

> 25 July 2024

- Add citation [`c1c3264`](https://github.com/emit-sds/emit-sds-l1b/commit/c1c3264bd127c48c02643a61fc297f72def28c3e)

#### [v1.4.1](https://github.com/emit-sds/emit-sds-l1b/compare/v1.4.0...v1.4.1)

> 1 September 2023

- Merge develop into main for v1.4.1 [`#20`](https://github.com/emit-sds/emit-sds-l1b/pull/20)
- Pass in rdn runconfig to get ffupdate paths [`9dfdd73`](https://github.com/emit-sds/emit-sds-l1b/commit/9dfdd7336c8eaca08d7c7643f99d206f3a2bb262)
- Update change log [`d135145`](https://github.com/emit-sds/emit-sds-l1b/commit/d1351450c16d59455864d779504475432ed0581b)

#### [v1.4.0](https://github.com/emit-sds/emit-sds-l1b/compare/v1.3.6...v1.4.0)

> 16 May 2023

- Merge develop into main for v1.4.0 [`#19`](https://github.com/emit-sds/emit-sds-l1b/pull/19)
- Wrapper script for new destriping [`#18`](https://github.com/emit-sds/emit-sds-l1b/pull/18)
- new flatfield utilities [`#17`](https://github.com/emit-sds/emit-sds-l1b/pull/17)
- Add calls to calibrate and destripe [`c78c720`](https://github.com/emit-sds/emit-sds-l1b/commit/c78c72049336fd20bf43543fbceb70edb524f885)
- Define paths in wrapper to keep runconfig cleaner. [`880caf1`](https://github.com/emit-sds/emit-sds-l1b/commit/880caf1422ea089a4d4cb84b23e610c72e0a2301)
- Only perform median and applyflat if we have more than 100 recent flat fields [`c1b348a`](https://github.com/emit-sds/emit-sds-l1b/commit/c1b348a4ed2f3ee0bb3603037883594db9d2ba27)

#### [v1.3.6](https://github.com/emit-sds/emit-sds-l1b/compare/v1.3.5...v1.3.6)

> 17 February 2023

- Merge develop into main for v1.3.6 [`#16`](https://github.com/emit-sds/emit-sds-l1b/pull/16)
- initial commit [`499a0d8`](https://github.com/emit-sds/emit-sds-l1b/commit/499a0d86b0d7883b4485b891227d66d458a13334)
- fixed edge effects [`e42abc5`](https://github.com/emit-sds/emit-sds-l1b/commit/e42abc5bb9c7cb035e0f89621a348b1bd87a2f95)
- Update change log [`3b5dc99`](https://github.com/emit-sds/emit-sds-l1b/commit/3b5dc99115a7d60e69a632251810d16f1916bf2c)

#### [v1.3.5](https://github.com/emit-sds/emit-sds-l1b/compare/v1.3.4...v1.3.5)

> 25 January 2023

- Merge develop into main for v1.3.5 [`#15`](https://github.com/emit-sds/emit-sds-l1b/pull/15)
- radiometric tweak [`#14`](https://github.com/emit-sds/emit-sds-l1b/pull/14)
- removed 680 change [`62cbea9`](https://github.com/emit-sds/emit-sds-l1b/commit/62cbea94b7aa2cd514264ea1e4233dc202fe0fb2)
- Update change log [`32234ef`](https://github.com/emit-sds/emit-sds-l1b/commit/32234ef1437cd777d1c9d833f330c41d7494fea6)

#### [v1.3.4](https://github.com/emit-sds/emit-sds-l1b/compare/v1.3.3...v1.3.4)

> 6 January 2023

- Merge develop into main for v1.3.4 [`#13`](https://github.com/emit-sds/emit-sds-l1b/pull/13)
- Recognize bad data flags [`#12`](https://github.com/emit-sds/emit-sds-l1b/pull/12)
- Fix NetCDF metadata issues [`5f2847a`](https://github.com/emit-sds/emit-sds-l1b/commit/5f2847ace4a888b62b7c6e403a66407a5183e7ff)
- Update change log [`b756425`](https://github.com/emit-sds/emit-sds-l1b/commit/b756425f4c44adfc4f58e7efc17af2e6664eac5c)

#### [v1.3.3](https://github.com/emit-sds/emit-sds-l1b/compare/v1.3.2...v1.3.3)

> 11 November 2022

- Merge develop into main for v1.3.3 [`#11`](https://github.com/emit-sds/emit-sds-l1b/pull/11)
- Update change log [`ea488c3`](https://github.com/emit-sds/emit-sds-l1b/commit/ea488c3cfe520037d4b8e6dd7bcc01d262c1ace2)
- add safeguards to flatfield [`de5c329`](https://github.com/emit-sds/emit-sds-l1b/commit/de5c329aaedf047ce78cf198bf3d28a6233472c1)

#### [v1.3.2](https://github.com/emit-sds/emit-sds-l1b/compare/v1.3.1...v1.3.2)

> 10 November 2022

- Merge develop into main for v1.3.2 [`#10`](https://github.com/emit-sds/emit-sds-l1b/pull/10)
- update fit flat field to handle saturation and bright pixels [`#9`](https://github.com/emit-sds/emit-sds-l1b/pull/9)
- Update change log [`fc6e8d8`](https://github.com/emit-sds/emit-sds-l1b/commit/fc6e8d8cb12e6f30b6a42f70051aa0529be770ad)
- bring hardcoded 4std into arguments [`c6e078e`](https://github.com/emit-sds/emit-sds-l1b/commit/c6e078eaa9ff48b68edd8f6c909b55fb719ef852)

#### [v1.3.1](https://github.com/emit-sds/emit-sds-l1b/compare/v1.3.0...v1.3.1)

> 7 November 2022

- Merge develop into main for v1.3.1 [`#8`](https://github.com/emit-sds/emit-sds-l1b/pull/8)
- multiband edge detection [`f67e5d8`](https://github.com/emit-sds/emit-sds-l1b/commit/f67e5d8937fc11550a35be0959177a48f852841f)
- Update change log [`564f464`](https://github.com/emit-sds/emit-sds-l1b/commit/564f464ad2a9ff6c5119734f9e5cb2210edb0d6d)

#### [v1.3.0](https://github.com/emit-sds/emit-sds-l1b/compare/v1.2.0...v1.3.0)

> 17 October 2022

- Merge develop into main for v1.3.0 [`#7`](https://github.com/emit-sds/emit-sds-l1b/pull/7)
- Merge IOC hotfixes into develop [`#6`](https://github.com/emit-sds/emit-sds-l1b/pull/6)
- Merge develop to hotfix v1.2.0 for deploy to ops [`#5`](https://github.com/emit-sds/emit-sds-l1b/pull/5)
- Merge develop to hotfix v1.2.0 for deploy to ops [`#4`](https://github.com/emit-sds/emit-sds-l1b/pull/4)
- Merging develop branch to hotfix 1.2.0 branch to deploy to ops [`#3`](https://github.com/emit-sds/emit-sds-l1b/pull/3)
- Updated standard deviation calculation [`#2`](https://github.com/emit-sds/emit-sds-l1b/pull/2)
- cubic interpolation for OSF Seams courtesy Phil B. [`0f11eca`](https://github.com/emit-sds/emit-sds-l1b/commit/0f11eca4b6416004248acf86b2cf963edf44c3a2)
- radiometry fix [`7ea80b6`](https://github.com/emit-sds/emit-sds-l1b/commit/7ea80b6948190a5ac085417613209a23c9ba6335)
- update to the blue end [`4028d2e`](https://github.com/emit-sds/emit-sds-l1b/commit/4028d2e5cfcfa36f02bde72a815bbb1283071325)

#### [v1.2.0](https://github.com/emit-sds/emit-sds-l1b/compare/v1.1.0...v1.2.0)

> 6 June 2022

- Merge develop into main for v1.2.0 [`#1`](https://github.com/emit-sds/emit-sds-l1b/pull/1)
- updated wavelengths and unit test [`aab2d49`](https://github.com/emit-sds/emit-sds-l1b/commit/aab2d49ca11acddaabf828bcdd485517b0feaeef)
- yet more cleanup [`30e57d9`](https://github.com/emit-sds/emit-sds-l1b/commit/30e57d981d8a59ddc247258a4f6be2e364db175e)
- final preflight calibrations? [`48ea986`](https://github.com/emit-sds/emit-sds-l1b/commit/48ea98697dfeeb9c38aa6906b9d402990e0a34bb)

#### [v1.1.0](https://github.com/emit-sds/emit-sds-l1b/compare/v1.0.0...v1.1.0)

> 4 May 2022

- Initial release (#4) [`#7`](https://github.com/emit-sds/emit-sds-l1b/pull/7)
- Removing predevelop branch [`#5`](https://github.com/emit-sds/emit-sds-l1b/pull/5)
- Netcdf [`#3`](https://github.com/emit-sds/emit-sds-l1b/pull/3)
- Merge pull request #6 from emit-sds/develop [`2250fb7`](https://github.com/emit-sds/emit-sds-l1b/commit/2250fb70cf9687efcb6e304d6c7c68394d8b5310)
- pointwise ghosts [`edd4f7d`](https://github.com/emit-sds/emit-sds-l1b/commit/edd4f7de78d24e1b964eeecad8f48cd60fee336b)
- new wavelength data [`82bf542`](https://github.com/emit-sds/emit-sds-l1b/commit/82bf542648aa1774da28b49a57bf57520d352d99)

### [v1.0.0](https://github.com/emit-sds/emit-sds-l1b/compare/v0.1.0...v1.0.0)

> 9 February 2022

- Initial release [`#4`](https://github.com/emit-sds/emit-sds-l1b/pull/4)

#### v0.1.0

> 26 January 2021

- Merge develop into main for Release 1 [`#2`](https://github.com/emit-sds/emit-sds-l1b/pull/2)
- Emit main integration [`#1`](https://github.com/emit-sds/emit-sds-l1b/pull/1)
- initial commit of unit test [`f61dab1`](https://github.com/emit-sds/emit-sds-l1b/commit/f61dab1d3522d7d9a5431006df869813cbefd339)
- initial checkin of synthetic cal files [`3758f9c`](https://github.com/emit-sds/emit-sds-l1b/commit/3758f9c83b7180fae53d13c8a9568e61c1fdba61)
- updates, bugfixes [`ae3be6a`](https://github.com/emit-sds/emit-sds-l1b/commit/ae3be6a6688f65b18a4340084475963113cdff6f)
