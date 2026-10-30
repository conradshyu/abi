# ABI Trace File Export

Written by [Conrad Shyu](mailto:conradshyu@hotmail.com)

> **Note:** This repository contains source code only. Example ABI trace files are not included.

## Overview

ABI trace files (also called chromatograms) use the ABIF tagged-data format. An ABIF file begins with a 128-byte header and stores its data in records identified by four-character tags. The source code in this repository parses those records and exports signal and peak data as CSV.

The tag descriptions below cover commonly encountered records; the list is not exhaustive.

| Tag | Data type | Description |
| --- | --- | --- |
| `DATA` | `short` | Raw or processed signal data for a channel |
| `Dye#` | `short` | Number of fluorescent dyes |
| `DyeN` | `string` | Fluorescent dye name |
| `DySN` | `string` | Fluorescent dye set name |
| `EVNT` | `string` | Instrument event, such as electrophoresis start or stop |
| `InSc` | `long` | Electrophoresis injection time in seconds |
| `InVt` | `long` | Electrophoresis injection voltage in volts |
| `LANE` | `short` | Capillary lane number |
| `MTRX` | `long` | Fluorescent-signal filter matrix |
| `PEAK` | `custom` | Fluorescent peak records |
| `PK_#` | `short` | Number of fluorescent peaks for a channel |
| `RUNT` | `time` | Electrophoresis start and stop times |
| `Sd#P` | `short` | Number of peaks from the size standard |
| `SpNm` | `string` | User-assigned sample name |
| `Std#` | `short` | Number of peaks defined by the size standard |
| `StdF` | `string` | Size-standard file name |
| `User` | `string` | Registered instrument user name |

## Build

From the repository directory, compile the command-line exporter with a C++ compiler:

```sh
g++ -I. abi2csv.cpp abifile.cpp abitag.cpp -o abi2csv
```

## Usage

Pass the file extension without a leading dot. For example, to process `.abi` files:

```sh
./abi2csv abi
```

The program searches for matching files in the current directory. On non-Windows builds it also searches subdirectories.

For each input file, the program writes:

- `<input>_raw.csv` — available CCD and GeneScan signal channels, with one column per channel.
- `<input>_peak.csv` — peak records for each available channel, with position, height, begin/end positions and heights, area, and size.

## Source files

| File | Purpose |
| --- | --- |
| `abi2csv.cpp` | Command-line driver and CSV writers |
| `abifile.cpp` | ABIF file loading, tag lookup, and signal/peak extraction |
| `abifile.h` | `AbiFile` class and signal/peak data structures |
| `abitag.cpp` | ABIF tag record parsing and data-type names |
| `abitag.h` | `AbiTagRecord` class |
| `.gitignore` | Ignores common compiler outputs |
| `LICENSE` | GNU General Public License, version 3 |
