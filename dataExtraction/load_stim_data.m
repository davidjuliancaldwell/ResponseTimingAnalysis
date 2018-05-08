function [stim,sing,tact,fsStim,fsSing,fsData,fsTact] = load_stim_data(Stim,Sing,ECO1,Tact)

stim = Stim.data;
fsStim = Stim.info.SamplingRateHz;

fsData = ECO1.info.SamplingRateHz;

sing = Sing.data;
fsSing = Sing.info.SamplingRateHz;

tact = Tact.data;
fsTact = Tact.info.SamplingRateHz;

end