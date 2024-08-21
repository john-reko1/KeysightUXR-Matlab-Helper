% Script here to test the functionality of the scope_interface script
%% Setup parameters
ip_address = '';  % Paste scope IP address
data_format = 'WORD';  % Only WORD and ASCII are supported
%% Set scope settings
scope = scope_interface;
% Note only one VISA device can be made on a single IP address
% If you are trying to set up another scope_interface, use different IP
% NEED to run initialize before any other functions
scope.initialize(ip_address) % Connect matlab to visa
scope.set_data_format(data_format) % set up either ASCII or WORD

%% Get X/Y axis data from scope
% Need to set x axis on scope, pull here for plotting
test_source = 'CHAN1'; % Use channel 1 as input source
[xInc, xOrg, xRef] = scope.get_x_axis(test_source); % x increment, origin, reference
[yInc, yOrg, yRef] = scope.get_y_axis(test_source); % x increment, origin, reference
xunits = scope.get_xunits(test_source); % x axis units (useful for debug)
yunits = scope.get_yunits(test_source); %x axis units (useful for debug)
%% Data type notes
% Support for two types of data formats as of 8/5/2024. If the type is
% WORD, then the data will be returned as double values depending on what
% analysis you are using (i.e. spectral). If data type is ASCII, then the
% output is a string with comma separated double values. This means that
% the user will need to use str2double(strsplit(data, ',')) to get actual
% data
%% Set/Get waveform analysis type
% Set waveform analysis type (i.e. spectral, nrz, pam, unspecified)
source1 = 'CHAN1'; % can be any function, channel, or waveform memory
analysis_type1 = 'SPECTRAL'; % Needs to be all caps
scope.set_analysis_type(source1, analysis_type1) % sets analysis type for source, changes on display

% For debugging, can also query analysis type from a source
type = scope.get_analysis_type(source); % Gets type in string format

%% Get generic waveform data (real valued)
source2 = 'CHAN1'; % Needs to be displayed on the scope
% Gets waveform from scope with set analysis type and returns its value
% If data format is WORD, returned is the numerical data
% If data format is ASCII, returned is a string

% WORD data format is simply like this
waveform = scope.get_waveform(source2); 

% ASCII format would look like this
scope.set_data_format('ASCII')
tmp_waveform_asc = scope.get_waveform(source2);
% string is 'val1, val2, ...', so split string by commas and convert
waveform_asc = str2double(strsplit(tmp_waveform_asc, ','));
%% Get IQ data
% Similar to getting arbitrary wave form, same concept
scope.set_data_format('WORD')
source3 = 'CHAN1';
iq_data = scope.get_iq(source3); % Sets to spectral analysis

% Alternatively, one can do the following
scope.set_analysis_type('SPECTRAL') % Also changes the display
raw_iq_data = scope.get_waveform(source3);
% Data is split "real, imaginary, real, imaginary", so split up
iq_data1 = raw_iq_data(1:2:end-1) + 1j * raw_iq_data(2:2:end); % Make sure sizes are equal
%% Function Notes
% Note that for all functions the format with the inputs is the same as
% defined in the FFT magnitude function section. The first argument is the
% source(s) (i.e. channels, wave memories, other functions). The second
% input is the function number. This is used when there are multiple
% functions that you want to document. If this input is left blank, then it
% defaults to function 1. The final input is used to determine if the
% function is displayed on the scope. If the function is not displayed,
% then measurements can be pulled to MATLAB. This variable should be set to
% 1 if it is the first time using a function as it will pause for 1 second
% so it can be displayed. After that, unless changing the function number,
% this variable should either be left out or set to another number besides
% one. Overall, the way to use these functions is as follows
% 
% x = scope.get_function(source, function_number, 1) % first time using
% x1 = scope.get_function(source, function_number) % 2nd time using
% x2 = scope.get_function(source, new_function_number, 1) % new function
%
% If function1 is already displaying, overwrite with
% x3 = scope.get_function(source) % defaults to function1
% 
% Overall, if one is taking a lot of measurements from one type of
% fucntion (i.e. fft_magnitude), then it is best to set up the function on
% the scope (either with the interface or the scope_interface function) and
% then read without the third input to the scope.get_function() functions.
% If there is only one loop (i.e. taking 10 measurements), then the code
% can look as follows:
% for k = 1:n
%   measurement(k, :) = scope.get_function(source, fun_num, k)
% 
% This will delay the first measurement, ensure the function is displayed,
% then allow the user to continuously download.
%% Get/Create FFT magnitude function
fun_num = 1; % function number
is_disp = 1; % First time using function, need to set to one for the first call
source4 = 'CHAN1'; % Take the FFT of channel 1

fft_mag0 = scope.get_fft_mag(source4, fun_num, is_disp);

% subsequent calls can ignore the is_disp variable as the function should
% be displayed on the scope, always use all 3 arguments when accessing a
% new function number
fft_mag1 = scope.get_fft_mag(source4, fun_num);

% Can also store in different function numberf
fun_num2 = 2;
fft_mag2 = scope.get_fft_mag(source4, fun_num2);

% Defaults to function 1 if the function number input is not provided
fft_mag3 = scope.get_fft_mag(source4); % Stored in function 1
%% Get/Create FFT phase function
% Example, same setup as FFT magnitude
% Takes the fft phase of channel 2, stores it in function 3
source5 = 'CHAN2';
fun_num_phase = 3;
fft_phase0 = scope.get_fft_phase(source5, fun_num_phase, 1); % First read
fft_phase1 = scope.get_fft_phase(source5, fun_num_phase); % Overwrites
%% Set X axis for FFT function
% Used to set x-axis to logarithmic/linear
% Only use if the function is set to FFT
x_axis_type = 'LIN'; % linear x axis, can also set to LOG
function_num = 2;
scope.set_fft_xaxis(x_axis_type, function_num);

% Defaults to function 1 if second input is blank
scope.set_fft_xaxis(x_axis_type) % For function1 if fft
%% Get/Create adding function
% Add two sources together
s1 = 'CHAN1';
s2 = 'WMEMORY1'; % Waveform memory slot 1
add_fun_number = 1;
is_disp = 1;

% Add the values from channel1 to wmemory1 and stores in function 1
add_data = scope.add_functions(s1, s2, add_fun_number, is_disp);
%% Get/Create PSD function
% Get the power spectral density for a input source
psd_source = 'CHAN1';
psd_fun_num = 1;
psd_disp = 1;

% Get the PSD for channel 1 and store+display with function 1
psd = scope.get_PSD(psd_source, psd_fun_num, psd_disp);
%% Get/Create frequency envelope function
% Finds frequency envelope of an amplitude modulated signal via Hilbert 
% transform, returned signal is the square root of the 
% real and imaginary parts of the signal
env_source = 'CHAN1';
env_fun_num = 2;
env_disp = 1;

% Get frequency envelope for channel1, store+display in function 2
envelope = scope.get_frequency_envelope(env_source, env_fun_num, env_disp);
%% Create Functions 
% The scope.display_function(source) function applies a function to a
% specified source. The available functions are: 
possible_funcs = ["Freq_Envelope", "Average", "Differentiate", "FFT_mag", ...
 "FFT_phase", "Integrate", "Low_pass", "High_pass", "Power_spectrum", "Invert"];
% These correspond to finding the frequency envelope, windowed average
% (over an inputted number of measurements), derivative, FFT magnitude, FFT
% angle, integration, low pass filtering (with 3 dB magnitude specified as
% an input), high pass (3 dB magnitude input), PSD, and inversion
% display the function on scope as such:
disp_source = 'CHAN1';
fun_name = "Power_spectrum";
fun_num = 2;
% Display PSD of channel 1 as function 2 on oscilloscope
scope.display_function(disp_source, fun_name, fun_num)
% Now use scope.get_waveform() function to get the raw data from the
% created function
disp_data = scope.get_waveform('FUNC2');
% This method should be used if many measurements from one function need to
% be made. Note that nothing is returned if you try to get data from a
% function that isn't displaying on the scope.
%% Query function type data/Clear function display
% One is able to see what each function is doing and clear the functions
% from the display

fun_query_num = 1;
fun_query = scope.query_function(fun_query_num); % String with function description

clear_nums = [1, 2, 4];
scope.clear_function_display(clear_nums) % Clears functions 1,2,4 from scope display

%% Display Real-Time Eye, measure jitter/eye width
% To get real time eye, one needs to set up an appropriate transmit signal
% that will yield the eye diagram they are looking for. This functionality
% will turn on the color grade and return a few values from the eye. This
% only works on channel 1 waveforms

scope.display_eye(true) % Displays eye on scope

% Measure jitter with either peak-to-peak (PP) or RMS methodology
method = 'PP'; % Other option is RMS
jitter = scope.eye_measure_jitter(method);

% Measure Eye width, works for anything with color gradient turned on
width = scope.eye_measure_width();

scope.display_eye(false) % Turns off real-time eye, color gradient is off
%% Save source to waveform memory
mem_source = 'CHAN1';
mem_loc = 2;
scope.save2memory(mem_source, mem_loc); % Save channel1 to wmemory 2
%% Write/Read arbitrary VISA commands
% If one wants to write other commands not listed above, there are 3
% functions to write arbitrary commands or output the visadev resource

% Write custom line to scope
scope.write(':ANALYZE:SIGNAL:TYPE? CHAN1') % Find the signal type for channel 1
scope.read() % Need to read after querying, same as readline() visadev command

v = scope.get_visa(); % Now have direct connection to scope and can use visadev commands
