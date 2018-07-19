function [IZC,AP,SPCC,UTMC,ZN] = s83tblspc(par)

% function [IZC,AP,SPCC,UTMC,ZN] = s83tblspc(par)
% convert geodesic to state plane co-ordinates
% original BY E. CARLSON / subroutine
% 15.03.12 M.P.

% CREATE THE STATE PLANE COORDINATE TABLES
%ZN = struct(135,1);  Char-Structure
SPCC = zeros(135,6);
IZC = zeros(135,1);
%AP = struct(135,1);  Char-Structure
UTMC = zeros(60,1);

% LOAD THE TABLE OF SPC STATE ZONE CODES
      IZC(1)=101;
      IZC(2)=102;
      IZC(3)=5001;
      IZC(4)=5002;
      IZC(5)=5003;
      IZC(6)=5004;
      IZC(7)=5005;
      IZC(8)=5006;
      IZC(9)=5007;
      IZC(10)=5008;
      IZC(11)=5009;
      IZC(12)=5010;
      IZC(13)=201;
      IZC(14)=202;
      IZC(15)=203;
      IZC(16)=301;
      IZC(17)=302;
      IZC(18)=401;
      IZC(19)=402;
      IZC(20)=403;
      IZC(21)=404;
      IZC(22)=405;
      IZC(23)=406;
      IZC(24)=501;
      IZC(25)=502;
      IZC(26)=503;
      IZC(27)=600;
      IZC(28)=700;
      IZC(29)=901;
      IZC(30)=902;
      IZC(31)=903;
      IZC(32)=1001;
      IZC(33)=1002;
      IZC(34)=5101;
      IZC(35)=5102;
      IZC(36)=5103;
      IZC(37)=5104;
      IZC(38)=5105;
      IZC(39)=1101;
      IZC(40)=1102;
      IZC(41)=1103;
      IZC(42)=1201;
      IZC(43)=1202;
      IZC(44)=1301;
      IZC(45)=1302;
      IZC(46)=1401;
      IZC(47)=1402;
      IZC(48)=1501;
      IZC(49)=1502;
      IZC(50)=1601;
      IZC(51)=1602;
      IZC(52)=1701;
      IZC(53)=1702;
      IZC(54)=1703;
      IZC(55)=1801;
      IZC(56)=1802;
      IZC(57)=1900;
      IZC(58)=2001;
      IZC(59)=2002;
      IZC(60)=0;
      IZC(61)=0;
      IZC(62)=0;
      IZC(63)=2111;
      IZC(64)=2112;
      IZC(65)=2113;
      IZC(66)=2201;
      IZC(67)=2202;
      IZC(68)=2203;
      IZC(69)=2301;
      IZC(70)=2302;
      IZC(71)=2401;
      IZC(72)=2402;
      IZC(73)=2403;
      IZC(74)=2500;
      IZC(75)=0;
      IZC(76)=0;
      IZC(77)=2600;
      IZC(78)=0;
      IZC(79)=2701;
      IZC(80)=2702;
      IZC(81)=2703;
      IZC(82)=2800;
      IZC(83)=2900;
      IZC(84)=3001;
      IZC(85)=3002;
      IZC(86)=3003;
      IZC(87)=3101;
      IZC(88)=3102;
      IZC(89)=3103;
      IZC(90)=3104;
      IZC(91)=3200;
      IZC(92)=3301;
      IZC(93)=3302;
      IZC(94)=3401;
      IZC(95)=3402;
      IZC(96)=3501;
      IZC(97)=3502;
      IZC(98)=3601;
      IZC(99)=3602;
      IZC(100)=3701;
      IZC(101)=3702;
      IZC(102)=3800;
      IZC(103)=3900;
      IZC(104)=4001;
      IZC(105)=4002;
      IZC(106)=4100;
      IZC(107)=4201;
      IZC(108)=4202;
      IZC(109)=4203;
      IZC(110)=4204;
      IZC(111)=4205;
      IZC(112)=4301;
      IZC(113)=4302;
      IZC(114)=4303;
      IZC(115)=4400;
      IZC(116)=4501;
      IZC(117)=4502;
      IZC(118)=4601;
      IZC(119)=4602;
      IZC(120)=4701;
      IZC(121)=4702;
      IZC(122)=4801;
      IZC(123)=4802;
      IZC(124)=4803;
      IZC(125)=4901;
      IZC(126)=4902;
      IZC(127)=4903;
      IZC(128)=4904;
      IZC(129)=5200;
      IZC(130)=0;
      IZC(131)=0;
      IZC(132)=5300;
      IZC(133)=5400;
% New Kentucky Single Zone
      IZC(134)=1600;
% New GUAM zone 5401 not on line yet
      IZC(135)=5401;

% LOAD THE PROPER TYPE OF PROJECTION
%          L=LAMBERT CONIC PROJECTION
%          T=TRANSVERSE MERCATOR PROJECTION
%          O=OBLIQUE MERCATOR PROJECTION
%          N=CONSTANTS NOT YET AVAILABLE (NEEDS PERIODIC UPDATING)

      AP(1).s='T';
      AP(2).s='T';
      AP(3).s='O';
      AP(4).s='T';
      AP(5).s='T';
      AP(6).s='T';
      AP(7).s='T';
      AP(8).s='T';
      AP(9).s='T';
      AP(10).s='T';
      AP(11).s='T';
      AP(12).s='L';
      AP(13).s='T';
      AP(14).s='T';
      AP(15).s='T';
      AP(16).s='L';
      AP(17).s='L';
      AP(18).s='L';
      AP(19).s='L';
      AP(20).s='L';
      AP(21).s='L';
      AP(22).s='L';
      AP(23).s='L';
      AP(24).s='L';
      AP(25).s='L';
      AP(26).s='L';
      AP(27).s='L';
      AP(28).s='T';
      AP(29).s='T';
      AP(30).s='T';
      AP(31).s='L';
      AP(32).s='T';
      AP(33).s='T';
      AP(34).s='T';
      AP(35).s='T';
      AP(36).s='T';
      AP(37).s='T';
      AP(38).s='T';
      AP(39).s='T';
      AP(40).s='T';
      AP(41).s='T';
      AP(42).s='T';
      AP(43).s='T';
      AP(44).s='T';
      AP(45).s='T';
      AP(46).s='L';
      AP(47).s='L';
      AP(48).s='L';
      AP(49).s='L';
      AP(50).s='L';
      AP(51).s='L';
      AP(52).s='L';
      AP(53).s='L';
      AP(54).s='L';
      AP(55).s='T';
      AP(56).s='T';
      AP(57).s='L';
      AP(58).s='L';
      AP(59).s='L';
      AP(60).s='N';
      AP(61).s='N';
      AP(62).s='N';
      AP(63).s='L';
      AP(64).s='L';
      AP(65).s='L';
      AP(66).s='L';
      AP(67).s='L';
      AP(68).s='L';
      AP(69).s='T';
      AP(70).s='T';
      AP(71).s='T';
      AP(72).s='T';
      AP(73).s='T';
      AP(74).s='L';
      AP(75).s='N';
      AP(76).s='N';
      AP(77).s='L';
      AP(78).s='N';
      AP(79).s='T';
      AP(80).s='T';
      AP(81).s='T';
      AP(82).s='T';
      AP(83).s='T';
      AP(84).s='T';
      AP(85).s='T';
      AP(86).s='T';
      AP(87).s='T';
      AP(88).s='T';
      AP(89).s='T';
      AP(90).s='L';
      AP(91).s='L';
      AP(92).s='L';
      AP(93).s='L';
      AP(94).s='L';
      AP(95).s='L';
      AP(96).s='L';
      AP(97).s='L';
      AP(98).s='L';
      AP(99).s='L';
      AP(100).s='L';
      AP(101).s='L';
      AP(102).s='T';
      AP(103).s='L';
      AP(104).s='L';
      AP(105).s='L';
      AP(106).s='L';
      AP(107).s='L';
      AP(108).s='L';
      AP(109).s='L';
      AP(110).s='L';
      AP(111).s='L';
      AP(112).s='L';
      AP(113).s='L';
      AP(114).s='L';
      AP(115).s='T';
      AP(116).s='L';
      AP(117).s='L';
      AP(118).s='L';
      AP(119).s='L';
      AP(120).s='L';
      AP(121).s='L';
      AP(122).s='L';
      AP(123).s='L';
      AP(124).s='L';
      AP(125).s='T';
      AP(126).s='T';
      AP(127).s='T';
      AP(128).s='T';
      AP(129).s='L';
      AP(130).s='N';
      AP(131).s='N';
      AP(132).s='N';
      AP(133).s='T';
      AP(134).s='L';
      AP(135).s='T';
      
% INITIALIZE CONSTANTS TABLES
for I=1:135
   for J=1:6
      SPCC(I,J)=0.00;
   end
end

% LOAD CONSTANTS BY EACH STATE APHABETICALLY
% TRANSVERSE MERCATOR WILL HAVE 4 CONSTANTS
%          1 - CENTRAL MERIDIAN (CM)
%          2 - FALSE EASTING VALUE AT THE CM (METERS)
%          3 - SOUTHERNMOST PARALLEL
%          4 - SCALE FACTOR
%          5 - FALSE NORTHING VALUE AT SOUTHERMOST PARALLEL (METERS)
% LAMBERT CONIC WILL HAVE 6 CONSTANTS
%          1 - C. M.
%          2 - FALSE EASTING AT CM (METERS)
%          3 - FALSE NORTHING FOR SOUTHERNMOST PARALLEL (METERS),
%              USUALLY EQUALS ZERO
%          4 - LATITUDE OF SO. STD. PARALLEL
%          5 - LATITUDE OF NO. STD. PARALLEL
%          6 - LATITUDE OF SOUTHERNMOST PARALLEL
% OBLIQUE MERCATOR HAS 6 CONSTANTS
%          1 - C. M.
%          2 - FALSE EASTING (METERS)
%          3 - FALSE NORTHING (METERS)
%          4 - AXIS AZIMUTH
%          5 - SOUTHERNMOST PARALLEL
%          6 - SCALE FACTOR

%              AL EAST
      SPCC(1,1)=s83tinline(85.00,50.00);
      SPCC(1,2)=200000.00;
      SPCC(1,3)=s83tinline(30.00,30.00);
      SPCC(1,4)=25000.00;
      SPCC(1,5)=0.00;
%              AL WEST
      SPCC(2,1)=s83tinline(87.00,30.00);
      SPCC(2,2)=600000.00;
      SPCC(2,3)=30.00;
      SPCC(2,4)=15000.00;
      SPCC(2,5)=0.00;
%              AK 1
      SPCC(3,1)=s83tinline(133.00,40.00);
      SPCC(3,2)=5000000.00;
      SPCC(3,3)=-5000000.00;
      SPCC(3,4)=atan(-0.7500);
      SPCC(3,5)=57.00;
      SPCC(3,6)=10000.00;
%              AK 2
      SPCC(4,1)=142.00;
      SPCC(4,2)=500000.00;
      SPCC(4,3)=54.00;
      SPCC(4,4)=10000.00;
      SPCC(4,5)=0.00;
%              AK 3
      SPCC(5,1)=146.00;
      SPCC(5,2)=500000.00;
      SPCC(5,3)=54.00;
      SPCC(5,4)=10000.00;
      SPCC(5,5)=0.00;
%              AK 4
      SPCC(6,1)=150.00;
      SPCC(6,2)=500000.00;
      SPCC(6,3)=54.00;
      SPCC(6,4)=10000.00;
      SPCC(6,5)=0.00;
%              AK 5
      SPCC(7,1)=154.00;
      SPCC(7,2)=500000.00;
      SPCC(7,3)=54.00;
      SPCC(7,4)=10000.00;
      SPCC(7,5)=0.00;
%              AK 6
      SPCC(8,1)=158.00;
      SPCC(8,2)=500000.00;
      SPCC(8,3)=54.00;
      SPCC(8,4)=10000.00;
      SPCC(8,5)=0.00;
%              AK 7
      SPCC(9,1)=162.00;
      SPCC(9,2)=500000.00;
      SPCC(9,3)=54.00;
      SPCC(9,4)=10000.00;
      SPCC(9,5)=0.00;
%              AK 8
      SPCC(10,1)=166.00;
      SPCC(10,2)=500000.00;
      SPCC(10,3)=54.00;
      SPCC(10,4)=10000.00;
      SPCC(10,5)=0.00;
%              AK 9
      SPCC(11,1)=170.00;
      SPCC(11,2)=500000.00;
      SPCC(11,3)=54.00;
      SPCC(11,4)=10000.00;
      SPCC(11,5)=0.00;
%              AK 10
      SPCC(12,1)=176.00;
      SPCC(12,2)=1000000.00;
      SPCC(12,3)=0.00;
      SPCC(12,4)=s83tinline(51.00,50.00);
      SPCC(12,5)=s83tinline(53.00,50.00);
      SPCC(12,6)=51.00;
%              AZ WEST
      SPCC(15,1)=s83tinline(113.00,45.00);
      SPCC(15,2)=213360.00;
      SPCC(15,3)=31.00;
      SPCC(15,4)=15000.00;
      SPCC(15,5)=0.00;
%              AZ CENTRAL
      SPCC(14,1)=s83tinline(111.00,55.00);
      SPCC(14,2)=213360.00;
      SPCC(14,3)=31.00;
      SPCC(14,4)=10000.00;
      SPCC(14,5)=0.00;
%              AZ EAST
      SPCC(13,1)=s83tinline(110.00,10.00);
      SPCC(13,2)=213360.00;
      SPCC(13,3)=31.00;
      SPCC(13,4)=10000.00;
      SPCC(13,5)=0.00;
%              AR NORTH
      SPCC(16,1)=92.00;
      SPCC(16,2)=400000.00;
      SPCC(16,3)=0.00;
      SPCC(16,4)=s83tinline(34.00,56.00);
      SPCC(16,5)=s83tinline(36.00,14.00);
      SPCC(16,6)=s83tinline(34.00,20.00);
%              AR SOUTH
      SPCC(17,1)=92.00;
      SPCC(17,2)=400000.00;
      SPCC(17,3)=400000.00;
      SPCC(17,4)=s83tinline(33.00,18.00);
      SPCC(17,5)=s83tinline(34.00,46.00);
      SPCC(17,6)=s83tinline(32.00,40.00);
%              CA 1
      SPCC(18,1)=122.00;
      SPCC(18,2)=2000000.00;
      SPCC(18,3)=500000.00;
      SPCC(18,4)=40.00;
      SPCC(18,5)=s83tinline(41.00,40.00);
      SPCC(18,6)=s83tinline(39.00,20.00);
%              CA 2
      SPCC(19,1)=122.00;
      SPCC(19,2)=2000000.00;
      SPCC(19,3)=500000.00;
      SPCC(19,4)=s83tinline(38.00,20.00);
      SPCC(19,5)=s83tinline(39.00,50.00);
      SPCC(19,6)=s83tinline(37.00,40.00);
      SPCC(20,1)=120.500;
%              CA 3
      SPCC(20,2)=2000000.00;
      SPCC(20,3)=500000.00;
      SPCC(20,4)=s83tinline(37.00,4.00);
      SPCC(20,5)=s83tinline(38.00,26.00);
      SPCC(20,6)=36.500;
%              CA 4
      SPCC(21,1)=119.00;
      SPCC(21,2)=2000000.00;
      SPCC(21,3)=500000.00;
      SPCC(21,4)=36.00;
      SPCC(21,5)=37.2500;
      SPCC(21,6)=s83tinline(35.00,20.00);
%              CA 5
      SPCC(22,1)=118.00;
      SPCC(22,2)=2000000.00;
      SPCC(22,3)=500000.00;
      SPCC(22,4)=s83tinline(34.00,2.00);
      SPCC(22,5)=s83tinline(35.00,28.00);
      SPCC(22,6)=33.500;
%              CA 6
      SPCC(23,1)=116.2500;
      SPCC(23,2)=2000000.00;
      SPCC(23,3)=500000.00;
      SPCC(23,4)=s83tinline(32.00,47.00);
      SPCC(23,5)=s83tinline(33.00,53.00);
      SPCC(23,6)=s83tinline(32.00,10.00);
      SPCC(24,1)=105.500;
%              CO NORTH
      SPCC(24,2)=914401.828900;
      SPCC(24,3)=304800.609600;
      SPCC(24,4)=s83tinline(39.00,43.00);
      SPCC(24,5)=s83tinline(40.00,47.00);
      SPCC(24,6)=s83tinline(39.00,20.00);
%              CO CENTRAL
      SPCC(25,1)=105.500;
      SPCC(25,2)=914401.828900;
      SPCC(25,3)=304800.609600;
      SPCC(25,4)=s83tinline(38.00,27.00);
      SPCC(25,5)=s83tinline(39.00,45.00);
      SPCC(25,6)=s83tinline(37.00,50.00);
%              CO SOUTH
      SPCC(26,1)=105.500;
      SPCC(26,2)=914401.828900;
      SPCC(26,3)=304800.609600;
      SPCC(26,4)=s83tinline(37.00,14.00);
      SPCC(26,5)=s83tinline(38.00,26.00);
      SPCC(26,6)=s83tinline(36.00,40.00);
      SPCC(27,1)=s83tinline(72.00,45.00);
%              CT
      SPCC(27,2)=304800.609600;
      SPCC(27,3)=152400.304800;
      SPCC(27,4)=s83tinline(41.00,12.00);
      SPCC(27,5)=s83tinline(41.00,52.00);
      SPCC(27,6)=s83tinline(40.00,50.00);
%              DE
      SPCC(28,1)=s83tinline(75.00,25.00);
      SPCC(28,2)=200000.00;
      SPCC(28,3)=38.00;
      SPCC(28,4)=200000.00;
      SPCC(28,5)=0.00;
%              FL EAST
      SPCC(29,1)=81.00;
      SPCC(29,2)=200000.00;
      SPCC(29,3)=s83tinline(24.00,20.00);
      SPCC(29,4)=17000.00;
      SPCC(29,5)=0.00;
%              FL WEST
      SPCC(30,1)=82.00;
      SPCC(30,2)=200000.00;
      SPCC(30,3)=s83tinline(24.00,20.00);
      SPCC(30,4)=17000.00;
      SPCC(30,5)=0.00;
%              FL NORTH
      SPCC(31,1)=s83tinline(84.00,30.00);
      SPCC(31,2)=600000.00;
      SPCC(31,3)=0.00;
      SPCC(31,4)=s83tinline(29.00,35.00);
      SPCC(31,5)=s83tinline(30.00,45.00);
      SPCC(31,6)=29.00;
      SPCC(32,1)=s83tinline(82.00,10.00);
%              GA EAST
      SPCC(32,2)=200000.00;
      SPCC(32,3)=30.00;
      SPCC(32,4)=10000.00;
      SPCC(32,5)=0.00;
%              GA WEST
      SPCC(33,1)=s83tinline(84.00,10.00);
      SPCC(33,2)=700000.00;
      SPCC(33,3)=30.00;
      SPCC(33,4)=10000.00;
      SPCC(33,5)=0.00;
%              HI 1
      SPCC(34,1)=s83tinline(155.00,30.00);
      SPCC(34,2)=500000.00;
      SPCC(34,3)=s83tinline(18.00,50.00);
      SPCC(34,4)=30000.00;
      SPCC(34,5)=0.00;
%              HI 2
      SPCC(35,1)=s83tinline(156.00,40.00);
      SPCC(35,2)=500000.00;
      SPCC(35,3)=s83tinline(20.00,20.00);
      SPCC(35,4)=30000.00;
      SPCC(35,5)=0.00;
%              HI 3
      SPCC(36,1)=158.00;
      SPCC(36,2)=500000.00;
      SPCC(36,3)=s83tinline(21.00,10.00);
      SPCC(36,4)=100000.00;
      SPCC(36,5)=0.00;
%              HI 4
      SPCC(37,1)=s83tinline(159.00,30.00);
      SPCC(37,2)=500000.00;
      SPCC(37,3)=s83tinline(21.00,50.00);
      SPCC(37,4)=100000.00;
      SPCC(37,5)=0.00;
%              HI 5
      SPCC(38,1)=s83tinline(160.00,10.00);
      SPCC(38,2)=500000.00;
      SPCC(38,3)=s83tinline(21.00,40.00);
      SPCC(38,4)=1.00;
      SPCC(38,5)=0.00;
%              ID EAST
      SPCC(39,1)=s83tinline(112.00,10.00);
      SPCC(39,2)=200000.00;
      SPCC(39,3)=s83tinline(41.00,40.00);
      SPCC(39,4)=19000.00;
      SPCC(39,5)=0.00;
%              ID CENTRAL
      SPCC(40,1)=114.00;
      SPCC(40,2)=500000.00;
      SPCC(40,3)=s83tinline(41.00,40.00);
      SPCC(40,4)=19000.00;
      SPCC(40,5)=0.00;
%              ID WEST
      SPCC(41,1)=s83tinline(115.00,45.00);
      SPCC(41,2)=800000.00;
      SPCC(41,3)=s83tinline(41.00,40.00);
      SPCC(41,4)=15000.00;
      SPCC(41,5)=0.00;
%              IL EAST
      SPCC(42,1)=s83tinline(88.00,20.00);
      SPCC(42,2)=300000.00;
      SPCC(42,3)=s83tinline(36.00,40.00);
      SPCC(42,4)=40000.00;
      SPCC(42,5)=0.00;
%              IL WEST
      SPCC(43,1)=s83tinline(90.00,10.00);
      SPCC(43,2)=700000.00;
      SPCC(43,3)=s83tinline(36.00,40.00);
      SPCC(43,4)=17000.00;
      SPCC(43,5)=0.00;
%              IN EAST
      SPCC(44,1)=s83tinline(85.00,40.00);
      SPCC(44,2)=100000.00;
      SPCC(44,3)=37.500;
      SPCC(44,4)=30000.00;
      SPCC(44,5)=250000.00;
%              IN WEST
      SPCC(45,1)=s83tinline(87.00,5.00);
      SPCC(45,2)=900000.00;
      SPCC(45,3)=37.500;
      SPCC(45,4)=30000.00;
      SPCC(45,5)=250000.00;
%              IA NORTH
      SPCC(46,1)=93.500;
      SPCC(46,2)=1500000.00;
      SPCC(46,3)=1000000.00;
      SPCC(46,4)=s83tinline(42.00,4.00);
      SPCC(46,5)=s83tinline(43.00,16.00);
      SPCC(46,6)=41.500;
%              IA SOUTH
      SPCC(47,1)=93.500;
      SPCC(47,2)=500000.00;
      SPCC(47,3)=0.00;
      SPCC(47,4)=s83tinline(40.00,37.00);
      SPCC(47,5)=s83tinline(41.00,47.00);
      SPCC(47,6)=40.00;
%              KS NORTH
      SPCC(48,1)=98.00;
      SPCC(48,2)=400000.00;
      SPCC(48,3)=0.00;
      SPCC(48,4)=s83tinline(38.00,43.00);
      SPCC(48,5)=s83tinline(39.00,47.00);
      SPCC(48,6)=s83tinline(38.00,20.00);
%              KS SOUTH
      SPCC(49,1)=98.500;
      SPCC(49,2)=400000.00;
      SPCC(49,3)=400000.00;
      SPCC(49,4)=s83tinline(37.00,16.00);
      SPCC(49,5)=s83tinline(38.00,34.00);
      SPCC(49,6)=s83tinline(36.00,40.00);
%              KY NORTH
      SPCC(50,1)=s83tinline(84.00,15.00);
      SPCC(50,2)=500000.00;
      SPCC(50,3)=0.00;
      SPCC(50,4)=s83tinline(37.00,58.00);
      SPCC(50,5)=s83tinline(38.00,58.00);
      SPCC(50,6)=37.500;
%              KY SOUTH
      SPCC(51,1)=s83tinline(85.00,45.00);
      SPCC(51,2)=500000.00;
      SPCC(51,3)=500000.00;
      SPCC(51,4)=s83tinline(36.00,44.00);
      SPCC(51,5)=s83tinline(37.00,56.00);
      SPCC(51,6)=s83tinline(36.00,20.00);
%              LA NORTH
      SPCC(52,1)=92.500;
      SPCC(52,2)=1000000.00;
      SPCC(52,3)=0.00;
      SPCC(52,4)=s83tinline(31.00,10.00);
      SPCC(52,5)=s83tinline(32.00,40.00);
      SPCC(52,6)=30.500;
%              LA S
      SPCC(53,1)=s83tinline(91.00,20.00);
      SPCC(53,2)=1000000.00;
      SPCC(53,3)=0.00;
      SPCC(53,4)=s83tinline(29.00,18.00);
      SPCC(53,5)=s83tinline(30.00,42.00);
      SPCC(53,6)=28.500;
%              LA OFF
      SPCC(54,1)=s83tinline(91.00,20.00);
      SPCC(54,2)=1000000.00;
      SPCC(54,3)=0.00;
      SPCC(54,4)=s83tinline(26.00,10.00);
      SPCC(54,5)=s83tinline(27.00,50.00);
      SPCC(54,6)=25.500;
%              ME EAST
      SPCC(55,1)=68.500;
      SPCC(55,2)=300000.00;
      SPCC(55,3)=s83tinline(43.00,40.00);
      SPCC(55,4)=10000.00;
      SPCC(55,5)=0.00;
%              ME WEST
      SPCC(56,1)=s83tinline(70.00,10.00);
      SPCC(56,2)=900000.00;
      SPCC(56,3)=s83tinline(42.00,50.00);
      SPCC(56,4)=30000.00;
      SPCC(56,5)=0.00;
%              MD
      SPCC(57,1)=77.00;
      SPCC(57,2)=400000.00;
      SPCC(57,3)=0.00;
      SPCC(57,4)=s83tinline(38.00,18.00);
      SPCC(57,5)=s83tinline(39.00,27.00);
      SPCC(57,6)=s83tinline(37.00,40.00);
%              MA M
      SPCC(58,1)=71.500;
      SPCC(58,2)=200000.00;
      SPCC(58,3)=750000.00;
      SPCC(58,4)=s83tinline(41.00,43.00);
      SPCC(58,5)=s83tinline(42.00,41.00);
      SPCC(58,6)=41.00;
%              MA ISLAND
      SPCC(59,1)=70.500;
      SPCC(59,2)=500000.00;
      SPCC(59,3)=0.00;
      SPCC(59,4)=s83tinline(41.00,17.00);
      SPCC(59,5)=s83tinline(41.00,29.00);
      SPCC(59,6)=41.00;
%              MI NORTH
      SPCC(63,1)=87.00;
      SPCC(63,2)=8000000.00;
      SPCC(63,3)=0.00;
      SPCC(63,4)=s83tinline(45.00,29.00);
      SPCC(63,5)=s83tinline(47.00,5.00);
      SPCC(63,6)=s83tinline(44.00,47.00);
%              MI CENTRAL
      SPCC(64,1)=s83tinline(84.00,22.00);
      SPCC(64,2)=6000000.00;
      SPCC(64,3)=0.00;
      SPCC(64,4)=s83tinline(44.00,11.00);
      SPCC(64,5)=s83tinline(45.00,42.00);
      SPCC(64,6)=s83tinline(43.00,19.00);
%              MI SOUTH
      SPCC(65,1)=s83tinline(84.00,22.00);
      SPCC(65,2)=4000000.00;
      SPCC(65,3)=0.00;
      SPCC(65,4)=s83tinline(42.00,06.00);
      SPCC(65,5)=s83tinline(43.00,40.00);
      SPCC(65,6)=41.500;
%              MN NORTH
      SPCC(66,1)=s83tinline(93.00,6.00);
      SPCC(66,2)=800000.00;
      SPCC(66,3)=100000.00;
      SPCC(66,4)=s83tinline(47.00,2.00);
      SPCC(66,5)=s83tinline(48.00,38.00);
      SPCC(66,6)=46.500;
%              MN CENTRAL
      SPCC(67,1)=s83tinline(94.00,15.00);
      SPCC(67,2)=800000.00;
      SPCC(67,3)=100000.00;
      SPCC(67,4)=s83tinline(45.00,37.00);
      SPCC(67,5)=s83tinline(47.00,3.00);
      SPCC(67,6)=45.00;
%              MN SOUTH
      SPCC(68,1)=94.00;
      SPCC(68,2)=800000.00;
      SPCC(68,3)=100000.00;
      SPCC(68,4)=s83tinline(43.00,47.00);
      SPCC(68,5)=s83tinline(45.00,13.00);
      SPCC(68,6)=43.00;
%              MS EAST
      SPCC(69,1)=s83tinline(88.00,50.00);
      SPCC(69,2)=300000.00;
      SPCC(69,3)=29.500;
      SPCC(69,4)=20000.00;
      SPCC(69,5)=0.00;
%              MS WEST
      SPCC(70,1)=s83tinline(90.00,20.00);
      SPCC(70,2)=700000.00;
      SPCC(70,3)=29.500;
      SPCC(70,4)=20000.00;
      SPCC(70,5)=0.00;
%              MO EAST
      SPCC(71,1)=90.500;
      SPCC(71,2)=250000.00;
      SPCC(71,3)=s83tinline(35.00,50.00);
      SPCC(71,4)=15000.00;
      SPCC(71,5)=0.00;
%              MO CENTRAL
      SPCC(72,1)=92.500;
      SPCC(72,2)=500000.00;
      SPCC(72,3)=s83tinline(35.00,50.00);
      SPCC(72,4)=15000.00;
      SPCC(72,5)=0.00;
%              MO WEST
      SPCC(73,1)=94.500;
      SPCC(73,2)=850000.00;
      SPCC(73,3)=s83tinline(36.00,10.00);
      SPCC(73,4)=17000.00;
      SPCC(73,5)=0.00;
%              MT
      SPCC(74,1)=s83tinline(109.00,30.00);
      SPCC(74,2)=600000.00;
      SPCC(74,3)=0.00;
      SPCC(74,4)=45.00;
      SPCC(74,5)=49.00;
      SPCC(74,6)=s83tinline(44.00,15.00);
%              NE
      SPCC(77,1)=100.00;
      SPCC(77,2)=500000.00;
      SPCC(77,3)=0.00;
      SPCC(77,4)=40.00;
      SPCC(77,5)=43.00;
      SPCC(77,6)=s83tinline(39.00,50.00);
%              NV EAST
      SPCC(79,1)=s83tinline(115.00,35.00);
      SPCC(79,2)=200000.00;
      SPCC(79,3)=34.7500;
      SPCC(79,4)=10000.00;
      SPCC(79,5)=8000000.00;
%              NV CENTRAL
      SPCC(80,1)=s83tinline(116.00,40.00);
      SPCC(80,2)=500000.00;
      SPCC(80,3)=34.7500;
      SPCC(80,4)=10000.00;
      SPCC(80,5)=6000000.00;
%              NV WEST
      SPCC(81,1)=s83tinline(118.00,35.00);
      SPCC(81,2)=800000.00;
      SPCC(81,3)=34.7500;
      SPCC(81,4)=10000.00;
      SPCC(81,5)=4000000.00;
%              NH
      SPCC(82,1)=s83tinline(71.00,40.00);
      SPCC(82,2)=300000.00;
      SPCC(82,3)=42.500;
      SPCC(82,4)=30000.00;
      SPCC(82,5)=0.00;
%              NJ
      SPCC(83,1)=74.500;
      SPCC(83,2)=150000.00;
      SPCC(83,3)=s83tinline(38.00,50.00);
      SPCC(83,4)=10000.00;
      SPCC(83,5)=0.00;
%              NM EAST
      SPCC(84,1)=s83tinline(104.00,20.00);
      SPCC(84,2)=165000.00;
      SPCC(84,3)=31.00;
      SPCC(84,4)=11000.00;
      SPCC(84,5)=0.00;
%              NM CENTRAL
      SPCC(85,1)=s83tinline(106.00,15.00);
      SPCC(85,2)=500000.00;
      SPCC(85,3)=31.00;
      SPCC(85,4)=10000.00;
      SPCC(85,5)=0.00;
%              NM WEST
      SPCC(86,1)=s83tinline(107.00,50.00);
      SPCC(86,2)=830000.00;
      SPCC(86,3)=31.00;
      SPCC(86,4)=12000.00;
      SPCC(86,5)=0.00;
%              NY EAST
      SPCC(87,1)=74.500;
      SPCC(87,2)=150000.00;
      SPCC(87,3)=s83tinline(38.00,50.00);
      SPCC(87,4)=10000.00;
      SPCC(87,5)=0.00;
%              NY CENTRAL
      SPCC(88,1)=s83tinline(76.00,35.00);
      SPCC(88,2)=250000.00;
      SPCC(88,3)=40.00;
      SPCC(88,4)=16000.00;
      SPCC(88,5)=0.00;
%              NY WEST
      SPCC(89,1)=s83tinline(78.00,35.00);
      SPCC(89,2)=350000.00;
      SPCC(89,3)=40.00;
      SPCC(89,4)=16000.00;
      SPCC(89,5)=0.00;
%              NY LI
      SPCC(90,1)=74.00;
      SPCC(90,2)=300000.00;
      SPCC(90,3)=0.00;
      SPCC(90,4)=s83tinline(40.00,40.00);
      SPCC(90,5)=s83tinline(41.00,2.00);
      SPCC(90,6)=s83tinline(40.00,10.00);
%              NC
      SPCC(91,1)=79.00;
      SPCC(91,2)=609601.2200;
      SPCC(91,3)=0.00;
      SPCC(91,4)=s83tinline(34.00,20.00);
      SPCC(91,5)=s83tinline(36.00,10.00);
      SPCC(91,6)=33.7500;
%              ND NORTH
      SPCC(92,1)=100.500;
      SPCC(92,2)=600000.00;
      SPCC(92,3)=0.00;
      SPCC(92,4)=s83tinline(47.00,26.00);
      SPCC(92,5)=s83tinline(48.00,44.00);
      SPCC(92,6)=47.00;
%              ND SOUTH
      SPCC(93,1)=100.500;
      SPCC(93,2)=600000.00;
      SPCC(93,3)=0.00;
      SPCC(93,4)=s83tinline(46.00,11.00);
      SPCC(93,5)=s83tinline(47.00,29.00);
      SPCC(93,6)=s83tinline(45.00,40.00);
%              OH NORTH
      SPCC(94,1)=82.500;
      SPCC(94,2)=600000.00;
      SPCC(94,3)=0.00;
      SPCC(94,4)=s83tinline(40.00,26.00);
      SPCC(94,5)=s83tinline(41.00,42.00);
      SPCC(94,6)=s83tinline(39.00,40.00);
%              OH SOUTH
      SPCC(95,1)=82.500;
      SPCC(95,2)=600000.00;
      SPCC(95,3)=0.00;
      SPCC(95,4)=s83tinline(38.00,44.00);
      SPCC(95,5)=s83tinline(40.00,2.00);
      SPCC(95,6)=38.00;
%              OK NORTH
      SPCC(96,1)=98.00;
      SPCC(96,2)=600000.00;
      SPCC(96,3)=0.00;
      SPCC(96,4)=s83tinline(35.00,34.00);
      SPCC(96,5)=s83tinline(36.00,46.00);
      SPCC(96,6)=35.00;
%              OK SOUTH
      SPCC(97,1)=98.00;
      SPCC(97,2)=600000.00;
      SPCC(97,3)=0.00;
      SPCC(97,4)=s83tinline(33.00,56.00);
      SPCC(97,5)=s83tinline(35.00,14.00);
      SPCC(97,6)=s83tinline(33.00,20.00);
%              OR NORTH
      SPCC(98,1)=120.500;
      SPCC(98,2)=2500000.00;
      SPCC(98,3)=0.00;
      SPCC(98,4)=s83tinline(44.00,20.00);
      SPCC(98,5)=46.00;
      SPCC(98,6)=s83tinline(43.00,40.00);
%              OR SOUTH
      SPCC(99,1)=120.500;
      SPCC(99,2)=1500000.00;
      SPCC(99,3)=0.00;
      SPCC(99,4)=s83tinline(42.00,20.00);
      SPCC(99,5)=44.00;
      SPCC(99,6)=s83tinline(41.00,40.00);
%              PA NORTH
      SPCC(100,1)=s83tinline(77.00,45.00);
      SPCC(100,2)=600000.00;
      SPCC(100,3)=0.00;
      SPCC(100,4)=s83tinline(40.00,53.00);
      SPCC(100,5)=s83tinline(41.00,57.00);
      SPCC(100,6)=s83tinline(40.00,10.00);
%              PA SOUTH
      SPCC(101,1)=s83tinline(77.00,45.00);
      SPCC(101,2)=600000.00;
      SPCC(101,3)=0.00;
      SPCC(101,4)=s83tinline(39.00,56.00);
      SPCC(101,5)=s83tinline(40.00,58.00);
      SPCC(101,6)=s83tinline(39.00,20.00);
%              RI
      SPCC(102,1)=71.500;
      SPCC(102,2)=100000.00;
      SPCC(102,3)=s83tinline(41.00,5.00);
      SPCC(102,4)=160000.00;
      SPCC(102,5)=0.00;
%              SC
      SPCC(103,1)=81.00;
      SPCC(103,2)=609600.00;
      SPCC(103,3)=0.00;
      SPCC(103,4)=32.500;
      SPCC(103,5)=s83tinline(34.00,50.00);
      SPCC(103,6)=s83tinline(31.00,50.00);
%              SD NORTH
      SPCC(104,1)=100.00;
      SPCC(104,2)=600000.00;
      SPCC(104,3)=0.00;
      SPCC(104,4)=s83tinline(44.00,25.00);
      SPCC(104,5)=s83tinline(45.00,41.00);
      SPCC(104,6)=s83tinline(43.00,50.00);
%              SD SOUTH
      SPCC(105,1)=s83tinline(100.00,20.00);
      SPCC(105,2)=600000.00;
      SPCC(105,3)=0.00;
      SPCC(105,4)=s83tinline(42.00,50.00);
      SPCC(105,5)=s83tinline(44.00,24.00);
      SPCC(105,6)=s83tinline(42.00,20.00);
%              TN
      SPCC(106,1)=86.00;
      SPCC(106,2)=600000.00;
      SPCC(106,3)=0.00;
      SPCC(106,4)=s83tinline(35.00,15.00);
      SPCC(106,5)=s83tinline(36.00,25.00);
      SPCC(106,6)=s83tinline(34.00,20.00);
%              TX NORTH
      SPCC(107,1)=101.500;
      SPCC(107,2)=200000.00;
      SPCC(107,3)=1000000.00;
      SPCC(107,4)=s83tinline(34.00,39.00);
      SPCC(107,5)=s83tinline(36.00,11.00);
      SPCC(107,6)=34.00;
%              TX NCENTRAL
      SPCC(108,1)=98.500;
      SPCC(108,2)=600000.00;
      SPCC(108,3)=2000000.00;
      SPCC(108,4)=s83tinline(32.00,8.00);
      SPCC(108,5)=s83tinline(33.00,58.00);
      SPCC(108,6)=s83tinline(31.00,40.00);
%              TX C
      SPCC(109,1)=s83tinline(100.00,20.00);
      SPCC(109,2)=700000.00;
      SPCC(109,3)=3000000.00;
      SPCC(109,4)=s83tinline(30.00,7.00);
      SPCC(109,5)=s83tinline(31.00,53.00);
      SPCC(109,6)=s83tinline(29.00,40.00);
%              TX SCENTRAL
      SPCC(110,1)=99.00;
      SPCC(110,2)=600000.00;
      SPCC(110,3)=4000000.00;
      SPCC(110,4)=s83tinline(28.00,23.00);
      SPCC(110,5)=s83tinline(30.00,17.00);
      SPCC(110,6)=s83tinline(27.00,50.00);
%              TX S
      SPCC(111,1)=98.500;
      SPCC(111,2)=300000.00;
      SPCC(111,3)=5000000.00;
      SPCC(111,4)=s83tinline(26.00,10.00);
      SPCC(111,5)=s83tinline(27.00,50.00);
      SPCC(111,6)=s83tinline(25.00,40.00);
%              UT NORTH
      SPCC(112,1)=111.500;
      SPCC(112,2)=500000.00;
      SPCC(112,3)=1000000.00;
      SPCC(112,4)=s83tinline(40.00,43.00);
      SPCC(112,5)=s83tinline(41.00,47.00);
      SPCC(112,6)=s83tinline(40.00,20.00);
%              UT CENTRAL
      SPCC(113,1)=111.500;
      SPCC(113,2)=500000.00;
      SPCC(113,3)=2000000.00;
      SPCC(113,4)=s83tinline(39.00,1.00);
      SPCC(113,5)=s83tinline(40.00,39.00);
      SPCC(113,6)=s83tinline(38.00,20.00);
%              UT SOUTH
      SPCC(114,1)=111.500;
      SPCC(114,2)=500000.00;
      SPCC(114,3)=3000000.00;
      SPCC(114,4)=s83tinline(37.00,13.00);
      SPCC(114,5)=s83tinline(38.00,21.00);
      SPCC(114,6)=s83tinline(36.00,40.00);
%              VT
      SPCC(115,1)=s83tinline(72.00,30.00);
      SPCC(115,2)=500000.00;
      SPCC(115,3)=s83tinline(42.00,30.00);
      SPCC(115,4)=28000.00;
      SPCC(115,5)=0.00;
%              VA NORTH
      SPCC(116,1)=78.500;
      SPCC(116,2)=3500000.00;
      SPCC(116,3)=2000000.00;
      SPCC(116,4)=s83tinline(38.00,2.00);
      SPCC(116,5)=s83tinline(39.00,12.00);
      SPCC(116,6)=s83tinline(37.00,40.00);
%              VA SOUTH
      SPCC(117,1)=78.500;
      SPCC(117,2)=3500000.00;
      SPCC(117,3)=1000000.00;
      SPCC(117,4)=s83tinline(36.00,46.00);
      SPCC(117,5)=s83tinline(37.00,58.00);
      SPCC(117,6)=s83tinline(36.00,20.00);
%              WA NORTH
      SPCC(118,1)=s83tinline(120.00,50.00);
      SPCC(118,2)=500000.00;
      SPCC(118,3)=0.00;
      SPCC(118,4)=47.500;
      SPCC(118,5)=s83tinline(48.00,44.00);
      SPCC(118,6)=47.00;
%              WA SOUTH
      SPCC(119,1)=120.500;
      SPCC(119,2)=500000.00;
      SPCC(119,3)=0.00;
      SPCC(119,4)=s83tinline(45.00,50.00);
      SPCC(119,5)=s83tinline(47.00,20.00);
      SPCC(119,6)=s83tinline(45.00,20.00);
%              WV NORTH
      SPCC(120,1)=79.500;
      SPCC(120,2)=600000.00;
      SPCC(120,3)=0.00;
      SPCC(120,4)=39.00;
      SPCC(120,5)=40.2500;
      SPCC(120,6)=38.500;
%              WV SOUTH
      SPCC(121,1)=81.00;
      SPCC(121,2)=600000.00;
      SPCC(121,3)=0.00;
      SPCC(121,4)=s83tinline(37.00,29.00);
      SPCC(121,5)=s83tinline(38.00,53.00);
      SPCC(121,6)=37.00;
%              WI NORTH
      SPCC(122,1)=90.00;
      SPCC(122,2)=600000.00;
      SPCC(122,3)=0.00;
      SPCC(122,4)=s83tinline(45.00,34.00);
      SPCC(122,5)=s83tinline(46.00,46.00);
      SPCC(122,6)=s83tinline(45.00,10.00);
%              WI CENTRAL
      SPCC(123,1)=90.00;
      SPCC(123,2)=600000.00;
      SPCC(123,3)=0.00;
      SPCC(123,4)=s83tinline(44.00,15.00);
      SPCC(123,5)=45.500;
      SPCC(123,6)=s83tinline(43.00,50.00);
%              WI SOUTH
      SPCC(124,1)=90.00;
      SPCC(124,2)=600000.00;
      SPCC(124,3)=0.00;
      SPCC(124,4)=s83tinline(42.00,44.00);
      SPCC(124,5)=s83tinline(44.00,4.00);
      SPCC(124,6)=42.00;
%              WY E
      SPCC(125,1)=s83tinline(105.00,10.00);
      SPCC(125,2)=200000.00;
      SPCC(125,3)=s83tinline(40.00,30.00);
      SPCC(125,4)=16000.00;
      SPCC(125,5)=0.00;
%              WY EC
      SPCC(126,1)=s83tinline(107.00,20.00);
      SPCC(126,2)=400000.00;
      SPCC(126,3)=s83tinline(40.00,30.00);
      SPCC(126,4)=16000.00;
      SPCC(126,5)=100000.00;
%              WY WC
      SPCC(127,1)=s83tinline(108.00,45.00);
      SPCC(127,2)=600000.00;
      SPCC(127,3)=s83tinline(40.00,30.00);
      SPCC(127,4)=16000.00;
      SPCC(127,5)=0.00;
%              WY W
      SPCC(128,1)=s83tinline(110.00,05.00);
      SPCC(128,2)=800000.00;
      SPCC(128,3)=s83tinline(40.00,30.00);
      SPCC(128,4)=16000.00;
      SPCC(128,5)=100000.00;

%             PUERTO RICO AND VIRGIN ISLANDS
      SPCC(129,1)=s83tinline(66.00,26.00);
      SPCC(129,2)=200000.00;
      SPCC(129,3)=200000.00;
      SPCC(129,4)=s83tinline(18.00,02.00);
      SPCC(129,5)=s83tinline(18.00,26.00);
      SPCC(129,6)=s83tinline(17.00,50.00);

%              GUAM
%     SPCC(133,1)=213.000;
%     SPCC(133,2)=500000.000;
%     SPCC(133,3)=0.000;
%     SPCC(133,4)=2500.000;
%     SPCC(133,5)=0.000;
%  
      SPCC(133,1)=s83tinline(215.00,15.00);
      SPCC(133,2)=100000.00;
      SPCC(133,3)=s83tinline(13.00,30.00);
      SPCC(133,4)=1.00;
      SPCC(133,5)=200000.00;

%              KY ONE
      SPCC(134,1)=85.7500;
      SPCC(134,2)=1500000.00;
      SPCC(134,3)=1000000.00;
      SPCC(134,4)=s83tinline(37.00,05.00);
      SPCC(134,5)=s83tinline(38.00,40.00);
      SPCC(134,6)=s83tinline(36.00,20.00);

%              GUAM New
      SPCC(135,1)=s83tinline(215.00,15.00);
      SPCC(135,2)=200000.00;
      SPCC(135,3)=s83tinline(13.00,30.00);
      SPCC(135,4)=1.00;
      SPCC(135,5)=100000.00;
      
%   UNIVERSAL TRANSVERSE MERCATOR HAS 4 CONSTANTS
%   LOAD CONSTANTS BY ZONES, 1 THRU 60
%            1 - CENTRAL MERIDIAN
%            2 - FALSE EASTING VALUE AT THE CM = 500,000.
%            3 - SOUTHERNMOST PARALLEL = 0.0
%            4 - SCALE FACTOR = 0.9996
%   SINCE THE LAST 3 CONSTANTS ARE ALWAYS THE SAME,
%         ONLY THE CENTRAL MERDIAN IS LOADED.

for I=1:60
   UTMC(I)=6.00*I-183.00;
end

      ZN(1).s='AL E';
      ZN(2).s='AL W';
      ZN(3).s='AK 1';
      ZN(4).s='AK 2';
      ZN(5).s='AK 3';
      ZN(6).s='AK 4';
      ZN(7).s='AK 5';
      ZN(8).s='AK 6';
      ZN(9).s='AK 7';
      ZN(10).s='AK 8';
      ZN(11).s='AK 9';
      ZN(12).s='AK10';
      ZN(13).s='AZ E';
      ZN(14).s='AZ C';
      ZN(15).s='AZ W';
      ZN(16).s='AR N';
      ZN(17).s='AR S';
      ZN(18).s='CA 1';
      ZN(19).s='CA 2';
      ZN(20).s='CA 3';
      ZN(21).s='CA 4';
      ZN(22).s='CA 5';
      ZN(23).s='CA 6';
      ZN(24).s='CO N';
      ZN(25).s='CO C';
      ZN(26).s='CO S';
      ZN(27).s='CT  ';
      ZN(28).s='DE  ';
      ZN(29).s='FL E';
      ZN(30).s='FL W';
      ZN(31).s='FL N';
      ZN(32).s='GA E';
      ZN(33).s='GA W';
      ZN(34).s='HI 1';
      ZN(35).s='HI 2';
      ZN(36).s='HI 3';
      ZN(37).s='HI 4';
      ZN(38).s='HI 5';
      ZN(39).s='ID E';
      ZN(40).s='ID C';
      ZN(41).s='ID W';
      ZN(42).s='IL E';
      ZN(43).s='IL W';
      ZN(44).s='IN E';
      ZN(45).s='IN W';
      ZN(46).s='IA N';
      ZN(47).s='IA S';
      ZN(48).s='KS N';
      ZN(49).s='KS S';
      ZN(50).s='KY N';
      ZN(51).s='KY S';
      ZN(52).s='LA N';
      ZN(53).s='LA S';
      ZN(54).s='LASH';
      ZN(55).s='ME E';
      ZN(56).s='ME W';
      ZN(57).s='MD  ';
      ZN(58).s='MA M';
      ZN(59).s='MA I';
      ZN(60).s='MI N';
      ZN(61).s='MI C';
      ZN(62).s='MI S';
      ZN(63).s='MI N';
      ZN(64).s='MI C';
      ZN(65).s='MI S';
      ZN(66).s='MN N';
      ZN(67).s='MN C';
      ZN(68).s='MN S';
      ZN(69).s='MS E';
      ZN(70).s='MS W';
      ZN(71).s='MO E';
      ZN(72).s='MO C';
      ZN(73).s='MO W';
      ZN(74).s='MT  ';
      ZN(75).s='MT  ';
      ZN(76).s='MT  ';
      ZN(77).s='NE  ';
      ZN(78).s='NE  ';
      ZN(79).s='NV E';
      ZN(80).s='NV C';
      ZN(81).s='NV W';
      ZN(82).s='NH  ';
      ZN(83).s='NJ  ';
      ZN(84).s='NM E';
      ZN(85).s='NM C';
      ZN(86).s='NM W';
      ZN(87).s='NY E';
      ZN(88).s='NY C';
      ZN(89).s='NY W';
      ZN(90).s='NY L';
      ZN(91).s='NC  ';
      ZN(92).s='ND N';
      ZN(93).s='ND S';
      ZN(94).s='OH N';
      ZN(95).s='OH S';
      ZN(96).s='OK N';
      ZN(97).s='OK S';
      ZN(98).s='OR N';
      ZN(99).s='OR S';
      ZN(100).s='PA N';
      ZN(101).s='PA S';
      ZN(102).s='RI  ';
      ZN(103).s='SC  ';
      ZN(104).s='SD N';
      ZN(105).s='SD S';
      ZN(106).s='TN  ';
      ZN(107).s='TX N';
      ZN(108).s='TXNC';
      ZN(109).s='TX C';
      ZN(110).s='TXSC';
      ZN(111).s='TX S';
      ZN(112).s='UT N';
      ZN(113).s='UT C';
      ZN(114).s='UT S';
      ZN(115).s='VT  ';
      ZN(116).s='VA N';
      ZN(117).s='VA S';
      ZN(118).s='WA N';
      ZN(119).s='WA S';
      ZN(120).s='WV N';
      ZN(121).s='WV S';
      ZN(122).s='WI N';
      ZN(123).s='WI C';
      ZN(124).s='WI S';
      ZN(125).s='WY E';
      ZN(126).s='WYEC';
      ZN(127).s='WYWC';
      ZN(128).s='WY W';
      ZN(129).s='PRVI';
      ZN(130).s='VIZ1';
      ZN(131).s='VISX';
      ZN(132).s='AS  ';
      ZN(133).s='GU  ';
      ZN(134).s='KY1Z';
      ZN(135).s='GU  ';

% disp('init done');
