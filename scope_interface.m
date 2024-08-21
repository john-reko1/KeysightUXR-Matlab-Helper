function scope = scope_interface
clear; clc; close all;
global initialized;
global format;
initialized = false;
format = 'ASCII';
scope.initialize = @initialize;
scope.set_data_format = @set_data_format;
scope.get_iq = @get_iq;
scope.get_fft_mag = @get_fft_mag;
scope.get_fft_phase = @get_fft_phase;
scope.set_fft_xaxis = @set_fft_xaxis;
scope.get_fft_peaks = @get_fft_peaks;
scope.get_xunits = @get_xunits;
scope.get_yunits = @get_yunits;
scope.get_x_axis = @get_x_axis;
scope.get_y_axis = @get_y_axis;
scope.display_eye = @display_eye;
scope.eye_measure_jitter = @eye_measure_jitter;
scope.eye_measure_width = @eye_measure_width;
scope.get_frequency_envelope = @get_frequency_envelope;
scope.query_function = @query_function;
scope.display_function = @create_function;
scope.add_functions = @add_functions;
scope.get_PSD = @get_PSD;
scope.get_analysis_type = @get_analysis_type;
scope.set_analysis_type = @set_analysis_type;
scope.clear_function_display = @clear_function_display;
scope.get_visa = @get_visa;
scope.save2memory = @save2memory;
scope.get_waveform = @get_waveform;
scope.write = @write_custom_command;
scope.read = @read;
scope.get_format = @get_format;
scope.close = @close_vsa;
end

function initialize(ip_address)
global initialized;
global v;

% Specify the VISA resource string of the oscilloscope
if nargin < 1
    ip_address = "TCPIP0::MJUXR0334B.lan::hislip0::INSTR" ;
end

% Create visadev object
v = visadev(ip_address);
v.Timeout = 10;

% Query the oscilloscope to ensure the connection is established
% Make sure to put check if the identificaiton is available
idn = writeread(v, '*IDN?'); % Find visa identification
disp(['Oscilloscope ID: ', idn]);
initialized = true;

% Clear Markers
writeline(v, ':MARKER1:ENABLE OFF');
writeline(v, ':MARKER1:SOURCE FUNCTION1');
writeline(v, ':MARKER1:TYPE RF');
writeline(v, ':MARKER2:ENABLE OFF');
writeline(v, ':MARKER2:SOURCE FUNCTION1');
writeline(v, ':MARKER2:TYPE RF');
writeline(v, ':MARKER3:ENABLE OFF');
writeline(v, ':MARKER3:SOURCE FUNCTION1');
writeline(v, ':MARKER3:TYPE RF');

% Re enable markers
writeline(v, ':MARKER1:ENABLE ON');
writeline(v, ':MARKER2:ENABLE ON');
writeline(v, ':MARKER3:ENABLE ON');
end

function set_data_format(dformat)
global initialized;
global v;
global format;

% All other functions are only set up for BYTE, WORD, and ASCII
% Can read other data formats if wanting to get the visadev object and
% write custom commands with functions write() and get_visa() functions

if ~initialized
    error('Need to run initialize function to connect VISA to scope')
end

possible_data = ["WORD", "BINARY", "BYTE", "ASCII", "FLOAT", "ASC", ...
    "BIN", "FLO"];

if nargin < 1
    valid_format = false; % If there are no input arguments
else
    tmp_format = upper(dformat);
    valid_format = ismember(tmp_format, possible_data);
end

if valid_format
    format = tmp_format;
    writeline(v, [':WAVEFORM:FORMAT ', tmp_format]);
    disp(['Data type is now ', tmp_format])
    if strcmp("WORD", tmp_format)
        writeline(v, ':WAVeform:BYTeorder LSBF')
    end
else
    disp('Not a valid data format type, default is ASCII')
end
end

function iq_data = get_iq(source)
global initialized;
global v;
global format;

if ~initialized
    error('Need to run initialize function to connect VISA to scope')
end

if nargin == 0
    source = 'CHAN1';
end

writeline(v, [':ANALyze:SIGNal:TYPE ', source,', SPECTRAL']);
% If not working, uncomment next two lines for debug
% writeline(v, [':ANALyze:SIGNal:TYPE? ', source])
% disp(readline(v))

% Set signal type to spectral for complex sampling
% Note that data still needs to be separated into real/imaginary
% components
writeline(v, [':WAVEFORM:SOURCE ', source])
if strcmp(format, "WORD")
    writeline(v, ':WAVEFORM:DATA?');
    iq_data = binblockread(v, 'short');
elseif strcmp(format, "BYTE")
    writeline(v, ':WAVEFORM:DATA?')
    iq_data = readbinblock(v, 'int8');
else
    iq_data = writeread(v, ':WAVEFORM:DATA?');
end
end

function xunits = get_xunits(source)
global v;
global initialized;
if ~initialized
    error('Need to run initialize function to connect VISA to scope')
end

if nargin == 0
    source = 'CHAN1';
end

writeline(v, [':WAVEFORM:SOURCE ', source])
xunits = strtrim(writeread(v, ':WAVEFORM:XUNITS?'));
end

function yunits = get_yunits(source)
global v;
global initialized;
if ~initialized
    error('Need to run initialize function to connect VISA to scope')
end

if nargin == 0
    source = 'CHAN1';
end

writeline(v, [':WAVEFORM:SOURCE ', source])
yunits = strtrim(writeread(v, ':WAVEFORM:YUNITS?'));
end

function [xInc, xOrg, xRef] = get_x_axis(source)
global v;
global initialized;

if ~initialized
    error('Need to run initialize function to connect VISA to scope')
end

if nargin == 0
    source = 'CHAN1';
end

writeline(v, [':WAVEFORM:SOURCE ', source])
writeline(v, ':WAVEFORM:XINCREMENT?');
xInc = str2double(readline(v));
writeline(v, ':WAVEFORM:XORIGIN?');
xOrg = str2double(readline(v));
writeline(v, ':WAVEFORM:XREFERENCE?');
xRef = str2double(readline(v));
end

function [yInc, yOrg, yRef] = get_y_axis(source)
global v;
global initialized;

if ~initialized
    error('Need to run initialize function to connect VISA to scope')
end

if nargin == 0
    source = 'CHAN1';
end

writeline(v, [':WAVEFORM:SOURCE ', source])
writeline(v, ':WAVEFORM:XINCREMENT?');
yInc = str2double(readline(v));
writeline(v, ':WAVEFORM:XORIGIN?');
yOrg = str2double(readline(v));
writeline(v, ':WAVEFORM:XREFERENCE?');
yRef = str2double(readline(v));
end

function fft_mag = get_fft_mag(source, tmp_fun_num, is_disp)
% Note that you have to turn the function on the display to work
global initialized;
global v;
global format;

if ~initialized
    error('Need to run initialize function to connect VISA to scope')
end
if nargin == 2
    is_disp = 0;
end

if nargin == 1
    tmp_fun_num = 1;
    is_disp = 0;
end

% Uses function 1 and 2 for mag and phase respectively
fun_num = num2str(tmp_fun_num);
writeline(v, [':FUNC', fun_num, ':FFTMagnitude ', source])
if is_disp == 1
    writeline(v, ['FUNC', fun_num, ':DISP 1'])
    pause(1)
end

% Query fft magnitude
writeline(v, [':WAVEFORM:SOURCE FUNC', fun_num])
if strcmp(format, "WORD")
    writeline(v, ':WAVEFORM:DATA?');
    fft_mag = binblockread(v, 'short'); % in general try new version
elseif strcmp(format, "BYTE")
    writeline(v, ':WAVEFORM:DATA?')
    fft_mag = readbinblock(v, 'int8');
else
    fft_mag = writeread(v, ':WAVEFORM:DATA?');
end
end

function fft_phase = get_fft_phase(source, tmp_fun_num, is_disp)
% Note that you have to turn the function on the display to work
global initialized;
global v;
global format;
if ~initialized
    error('Need to run initialize function to connect VISA to scope')
end

if nargin == 1
    tmp_fun_num = 1;
    is_disp = 0;
end

if nargin == 2
    is_disp = 0;
end

% Uses function 1 and 2 for mag and phase respectively
fun_num = num2str(tmp_fun_num);
writeline(v, [':FUNC', fun_num, ':FFTPhase ', source])
% Use for first time running function in program
if is_disp == 1
    writeline(v, ['FUNC', fun_num, ':DISP 1'])
    pause(1)
end
% Query fft magnitude
writeline(v, [':WAVEFORM:SOURCE FUNC', fun_num])
if strcmp(format, "WORD")
    writeline(v, ':WAVEFORM:DATA?');
    fft_phase = binblockread(v, 'short');
elseif strcmp(format, "BYTE")
    writeline(v, ':WAVEFORM:DATA?')
    fft_phase = readbinblock(v, 'int8');
else
    fft_phase = writeread(v, ':WAVEFORM:DATA?');
end
end

function set_fft_xaxis(linlog, fun_num)
global initialized;
global v;
if ~initialized
    error('Need to run initialize function to connect VISA to scope')
end
if nargin == 1
    fun_num = 1;
end
if ~strcmp(linlog, 'LOG') && ~strcmp(linlog, 'LIN')
    error('Needs to be LIN or LOG for linear or logarithmic')
else
    writeline(v, [':FUNCtion', num2str(fun_num), ':FFT:HSCale ' , linlog])
end
end
function [peak_mag, peak_freq] = get_fft_peaks(fun_num, num_peaks, is_disp)

global initialized;
global v;
if ~initialized
    error('Need to run initialize function to connect VISA to scope')
end
if nargin < 3
    is_disp = 0;
end
writeline(v, [':FUNC', num2str(fun_num), ':FFT:PEAK:COUNT ', num2str(num_peaks)])
% Turn on markers
writeline(v, [':FUNC', num2str(fun_num), ':FFT:PEAK:STATE ON'])
if is_disp == 1
    pause(1)
end
peak_mag = writeread(v, [':FUNC', num2str(fun_num), ':FFT:PEAK:MAGN?']);
peak_freq = writeread(v, [':FUNC', num2str(fun_num), ':FFT:PEAK:FREQ?']);
end
function data = add_functions(source1, source2, fun_num, is_disp)
global initialized;
global v;
global format;
if ~initialized
    error('Need to run initialize function to connect VISA to scope')
end
if nargin == 2
    fun_num = 1;
    is_disp = 0;
end

if nargin == 3
    is_disp = 0;
end
writeline(v, [':FUNC', num2str(fun_num), ':ADD ', source1, ', ' source2])
writeline(v, [':WAVEFORM:SOURCE FUNC', num2str(fun_num)])
if is_disp == 1
    writeline(v, [':FUNC', num2str(fun_num), ':DISP 1'])
    pause(1);
end
if strcmp(format, "WORD")
    writeline(v, ':WAVEFORM:DATA?');
    data = binblockread(v, 'short');
elseif strcmp(format, "BYTE")
    writeline(v, ':WAVEFORM:DATA?')
    data = readbinblock(v, 'int8');
else
    data = writeread(v, ':WAVEFORM:DATA?');
end
end
function psd = get_PSD(source, fun_num, is_disp)
global initialized;
global v;
global format;
if ~initialized
    error('Need to run initialize function to connect VISA to scope')
end
if nargin == 1
    fun_num = 1;
    is_disp = 0;
end
if nargin == 2
    is_disp = 0;
end
% set function to PSD
writeline(v, [':FUNC', num2str(fun_num), ':PSD ', source])
if is_disp == 1
    writeline(v, [':FUNC', num2str(fun_num), ':DISP 1'])
    pause(1);
end
writeline(v, [':WAVEFORM:SOURCE FUNC', num2str(fun_num)])
if strcmp(format, "WORD")
    writeline(v, ':WAVEFORM:DATA?');
    psd = binblockread(v, 'short');
elseif strcmp(format, "BYTE")
    writeline(v, ':WAVEFORM:DATA?')
    psd = readbinblock(v, 'int8');
else
    psd = writeread(v, ':WAVEFORM:DATA?');
end
end
function set_analysis_type(source, type)
global initialized;
global v;
if ~initialized
    error('Need to run initialize function to connect VISA to scope')
end
writeline(v, [':ANALyze:SIGNal:TYPE ', source,', ', type]);
end
function type = get_analysis_type(source)
global initialized;
global v;
if ~initialized
    error('Need to run initialize function to connect VISA to scope')
end
if nargin == 0
    source = 'CHAN1';
end
type = writeread(v, [':ANALyze:SIGNal:TYPE? ', source]);
end
function envelope = get_frequency_envelope(source, tmp_fun_num, is_disp)
global initialized;
global v;
global format;
if ~initialized
    error('Need to run initialize function to connect VISA to scope')
end

if nargin == 1
    tmp_fun_num = 1;
    is_disp = 0;
end
if nargin == 2
    is_disp = 0;
end
fun_num = num2str(tmp_fun_num);
writeline(v, [':FUNC', fun_num, ':ADEMod ', source])
if is_disp == 1
    writeline(v, [':FUNC', fun_num, ':DISPLAY 1'])
    pause(1);
end
writeline(v, [':WAVEFORM:SOURCE ', 'FUNC', fun_num])
if strcmp(format, "WORD")
    writeline(v, ':WAVEFORM:DATA?');
    envelope = binblockread(v, 'short');
elseif strcmp(format, "BYTE")
    writeline(v, ':WAVEFORM:DATA?')
    envelope = readbinblock(v, 'int8');
else
    envelope = writeread(v, ':WAVEFORM:DATA?');
end
end
function create_function(source, fun_name, fun_num)
global v;
global initialized;
if ~initialized
    error('Need to run initialize function to connect VISA to scope')
end
keys = ["Freq_Envelope", "Average", "Differentiate", "FFT_mag", ...
    "FFT_phase", "Integrate", "Low_pass", "High_pass", "Power_spectrum", "Invert"];
values = ["ADEM", "AVER", "DIFF", "FFTM", "FFTP", "INT", "LOWP", ...
    "HIGH", "PSD", "INV"];
valid_key = false;
for k = 1:size(keys, 2)
    if keys(k) == fun_name
        valid_key = true;
    end
end
if ~valid_key
    error('Function name is not allowed!!')
end
d = dictionary(keys, values);
if strcmp(fun_name, "Low_pass") || strcmp(fun_name, "High_pass")
    bw = input('Specify the 3dB bandwidth for the LPF: ');
    if ~isa(bw, 'double')
        error('input must be a double')
    end
    mag = floor(log(bw) / log(10));
    prefix = bw / 10^mag;
    tmp = convertStringsToChars(d(fun_name));
    writeline(v, [':FUNC', num2str(fun_num), ':', tmp, ' ', source, ...
        ',', num2str(prefix), 'E', num2str(mag)])
elseif strcmp(fun_name, "Average")
    % Averages 1:num_averages samples at first and every n reads after is
    % 1+n:num_averages+n samples, scope is default in powers of 2
    % can use for average number of samples
    num_averages = input('Number of samples for average: ');
    tmp = convertStringsToChars(d(fun_name));
    writeline(v, [':FUNC', num2str(fun_num), ':', tmp, ' ', source, ', ', num2str(num_averages)])
else
    tmp = convertStringsToChars(d(fun_name));
    writeline(v, [':FUNC', num2str(fun_num), ':', tmp, ' ', source])
end
writeline(v, [':FUNC', num2str(fun_num), ':DISP ON'])
end
function clear_function_display(fun_nums)
% Clears all functions with function number in input array
global initialized;
global v;
if ~initialized
    error('Need to run initialize function to connect VISA to scope')
end
for tmp_k = fun_nums
    k = num2str(tmp_k);
    %disp([':FUNC', k, ':DISP 0'])
    writeline(v, [':FUNC', k, ':DISP 0'])
end
end
function fun = query_function(function_num)
global initialized;
global v;
if ~initialized
    error('Need to run initialize function to connect VISA to scope')
end
% Turn on so you can see operator and operand
writeline(v, 'SYSTEM:HEADER 1')
fun = writeread(v, ['FUNC', num2str(function_num), '?']);
disp(fun)
end
function display_eye(on_off)
global v;
global initialized;
if ~initialized
    error('Need to run initialize function to connect VISA to scope')
end
%Activate real time eye for every signal with colorgrade activated
% look into beam focusing
if on_off
    writeline(v, 'MTEST:FOLDING 1')
else
    writeline(v, 'MTEST:FOLDING 0')
end
end
function jitter = eye_measure_jitter(format)
global v;
global initialized;
if ~initialized
    error('Need to run initialize function to connect VISA to scope')
end
writeline(v, [':MEASURE:CGRADE:JITTER ', format])
jitter = writeread(v, [':MEASURE:CGRADE:JITTER? ', format]);
end
function ewidth = eye_measure_width()
global v;
global initialized;
if ~initialized
    error('Need to run initialize function to connect VISA to scope')
end
writeline(v, ':MEASURE:CGRADE:EWIDTH')
ewidth = writeread(v, ':MEASURE:CGRADE:EWIDTH? ');
end
function save2memory(source, mem_loc)
global v;
global initialized;
if ~initialized
    error('Need to run initialize function to connect VISA to scope')
end
tmp_mem_loc = num2str(mem_loc);
writeline(v, [':WMEM', tmp_mem_loc, ':SAVE ', source])
end
function waveform = get_waveform(source)
global v;
global initialized;
global format;
if ~initialized
    error('Need to run initialize function to connect VISA to scope')
end
writeline(v, [':WAVEFORM:SOURCE ', source])
if strcmp(format, "WORD")
    writeline(v, ':WAVEFORM:DATA?');
    waveform = binblockread(v, 'short');
elseif strcmp(format, "BYTE")
    writeline(v, ':WAVEFORM:DATA?')
    waveform = readbinblock(v, 'int8');
else
    waveform = writeread(v, ':WAVEFORM:DATA?');
end
end
function v1 = get_visa()
global v;
global initialized;
if ~initialized
    error('Need to run initialize function to connect VISA to scope')
end
v1 = v;
end
function write_custom_command(command)
global v;
global initialized;
if ~initialized
    error('Need to run initialize function to connect VISA to scope')
end
writeline(v, command)
end
function q = read()
global v;
q = readline(v);
end
function form = get_format()
global format;
form = format;
end
function close_vsa()
global v;
fclose(v);
delete(v);
clear v;
end