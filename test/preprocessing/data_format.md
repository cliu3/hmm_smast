# Format of input data
The HMM Geolocation Toolbox supports raw ASCII data file from Star-ODDI DSTs downloaded using SeaStar.

The following is an example of a data file. Please configure SeaStar software properly so that the data file generated is as close to the required format as possible. For more information please refer to Chapter 8: Settings in the [Star-ODDI's User Manual for SeaStar](http://www.star-oddi.com/updates/SeaStar/DstTD.pdf). Specifically, we require:
* The first line of data to be in the 15th line of the file, and the first 14 lines to start with a number sign/pound sign `#`;
* Date and time format to be `dd.mm.yy HH:MM:SS`;
* Use decimal points (do not use decimal commas);
* Temperature values to be the third column of the data and in &deg;C;
* Depth values to be the fourth column of the data and in meters;
* At least two decimal places to be used if possible;
* Four underscores (`____`) to be used to indicate missing data.

        #0	Date-time:	5/28/2010 9:23:57 AM
        #1	Recorder:	1M11951
        #2	File type:	1
        #3	Columns:	4
        #4	Channels:	2
        #5	Field separation:	0
        #6	Decimal point:	1
        #7	Date def.:	0	0
        #8	Time def.:	0
        #9	Channel 1:	Temperature(°C)	Temp(°C)	2	1
        #10	Channel 2:	Depth(m)	Depth(m)	4	2
        #11	Reconvertion:	0
        #19	Line color:	1	2	3	4
        #30	Trend Type Number:	1
        1	07.05.10 05:00:00	17.30	0.8571
        2	07.05.10 05:15:00	____	0.9683
        3	07.05.10 05:30:00	____	1.4926
        4	07.05.10 05:45:00	____	1.9294
        5	07.05.10 06:00:00	____	2.3661
        6	07.05.10 06:15:00	____	2.7155
        7	07.05.10 06:30:00	____	2.8902
        8	07.05.10 06:45:00	____	3.2395
        9	07.05.10 07:00:00	____	25.9077
        10	07.05.10 07:15:00	____	59.4148
        11	07.05.10 07:30:00	____	59.5888
        12	07.05.10 07:45:00	4.58	54.4479
        13	07.05.10 08:00:00	____	54.5349
        14	07.05.10 08:15:00	____	54.4479
        15	07.05.10 08:30:00	____	54.0998
        16	07.05.10 08:45:00	____	54.3609
        17	07.05.10 09:00:00	____	54.3609
        18	07.05.10 09:15:00	____	54.2739
        19	07.05.10 09:30:00	____	54.2739
        20	07.05.10 09:45:00	____	54.1868
        21	07.05.10 10:00:00	____	54.0998
        22	07.05.10 10:15:00	____	54.0128
        23	07.05.10 10:30:00	4.54	53.8245
        24	07.05.10 10:45:00	____	53.5635
        25	07.05.10 11:00:00	____	53.4764
        26	07.05.10 11:15:00	____	53.3894
        27	07.05.10 11:30:00	____	53.2154
        28	07.05.10 11:45:00	____	53.2154
        29	07.05.10 12:00:00	____	53.4764
        30	07.05.10 12:15:00	____	53.3894
        31	07.05.10 12:30:00	____	53.3894
        32	07.05.10 12:45:00	____	53.3024
        33	07.05.10 13:00:00	____	53.3894
        34	07.05.10 13:15:00	4.58	53.3166
        35	07.05.10 13:30:00	____	53.4037
        36	07.05.10 13:45:00	____	53.4907
        37	07.05.10 14:00:00	____	53.4907
        38	07.05.10 14:15:00	____	53.5777
        39	07.05.10 14:30:00	____	53.7517
        40	07.05.10 14:45:00	____	53.8388
        41	07.05.10 15:00:00	____	53.9258
        42	07.05.10 15:15:00	____	54.0128
        43	07.05.10 15:30:00	____	54.0998
        44	07.05.10 15:45:00	____	54.1868
        45	07.05.10 16:00:00	4.58	54.2739
