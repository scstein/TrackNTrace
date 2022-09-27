function [sync, tcspc, chan, special, num, loc, head] = PTU_Read(name, cnts, head)
%
%  function [sync, tcspc, chan, special, num] = HT3_Read(name, cnts, head)
%
%  This function reads single-photon data from the file 'name'
%
%  If 'cnts' contains a number larger than 0, the routine reads 'cnts'
%  records the data stream or up the end of the file.
%
%  If 'cnts' contains two numbers [cnts(1) cnts(2)], the routine proceeds
%  to the position cnts(1) before readinf the cnts(2) records of data.
%
%  The output variables contain the followig data:
%  sync    : number of the sync events that preceeded this detection event
%  tcspc   : number of the tcspc-bin of the event
%  chan    : number of the input channel of the event (detector-number)
%  special : indicator of the event-type (0: photon; else : virtual photon)
%  num     : counter of the records that were actually read
%  loc     : number of overcounts after last valid photon

rtPicoHarpT3     = hex2dec('00010303');
rtPicoHarpT2     = hex2dec('00010203');
rtHydraHarpT3    = hex2dec('00010304');
rtHydraHarpT2    = hex2dec('00010204');
rtHydraHarp2T3   = hex2dec('01010304');
rtHydraHarp2T2   = hex2dec('01010204');
rtTimeHarp260NT3 = hex2dec('00010305');
rtTimeHarp260NT2 = hex2dec('00010205');
rtTimeHarp260PT3 = hex2dec('00010306');
rtTimeHarp260PT2 = hex2dec('00010206');
rtMultiHarpT3    = hex2dec('00010307');
rtMultiHarpT2    = hex2dec('00010207');


if (nargin<3)||(isempty(head))
    head = PTU_Read_Head(name);
end

if ~isempty(head)
    
    if (nargin<2)||isempty(cnts)
        cnts = [0 0];
    end;
    
    if numel(cnts)<2
        cnts = [0 cnts];
    end;

    
    if cnts(2)>0
        
        fid = fopen(name);
        
        if fid<1
            % Retry once
            pause(0.1);
            fid = fopen(name);
            warning('Could not open <%s>. On first try.\n', name);
        end 
        if fid<1
            warning('Could not open <%s>. Aborted.\n', name);
            return;
        else
                            
            fseek(fid, head.length, 'bof');
            
            if cnts(1)>1
                fsuc = fseek(fid, 4*(cnts(1)-1), 'cof');
                if fsuc<0
                    warning('Start read is greater than maximum number of photons.');
                    [sync, tcspc, chan, special] = deal([]);
                    num = 0;
                    loc = 0;
                    fclose(fid);
                    return;
                end                
            end
                                    
            [T3Record num] = fread(fid, cnts(2), 'ubit32'); % all 32 bits:
                        
            switch head.TTResultFormat_TTTRRecType;
                case rtPicoHarpT3
                    
                    WRAPAROUND=65536;
                    
                    sync    = bitand(T3Record,65535);              % the lowest 16 bits:
                    chan    = bitand(bitshift(T3Record,-28),15);   % the upper 4 bits:
                    tcspc   = bitand(bitshift(T3Record,-16),4095);                    
                    special = (chan==15).*bitand(tcspc,15);
                    
                    ind  = ((chan==15) & bitand(tcspc,15)==0);


                case rtPicoHarpT2

                    WRAPAROUND=210698240;
                    
                    sync    = bitand(T3Record,268435455);         %the lowest 28 bits
                    tcspc   = bitand(T3Record,15);                %the lowest 4 bits
                    chan    = bitand(bitshift(T3Record,-28),15);  %the next 4 bits
                    special = (chan==15).*bitand(tcspc,15);
                    
                    ind  = ((chan==15) & bitand(tcspc,15)==0);
                    
                case {rtHydraHarpT3, rtHydraHarp2T3, rtTimeHarp260NT3, rtTimeHarp260PT3, rtMultiHarpT3}
                    
                    WRAPAROUND=1024;

                    special = bitand(T3Record,2147483648)>0;
                    chan = bitand(bitshift(T3Record,-25),63);
                    tcspc = bitand(bitshift(T3Record,-10),32767);
                    sync = bitand(T3Record,1023);
                    
                    ind  = (special==1 & chan==63);

                    special = special.*chan;
                    
                case rtHydraHarpT2

                    WRAPAROUND=33552000;

                    sync    = bitand(T3Record,33554431);           % the last 25 bits
                    chan    = bitand(bitshift(T3Record,-25),63);   % the next 6 bits
                    tcspc   = bitand(chan,15);
                    special = bitand(bitshift(T3Record,-31),1);    % the last bit

                    ind  = (special==1 & chan==63);

                    special = special.*chan;
                                        
                case {rtHydraHarp2T2, rtTimeHarp260NT2, rtTimeHarp260PT2, rtMultiHarpT2}

                    WRAPAROUND=33554432;
                    
                    sync    = bitand(T3Record,33554431);           % the last 25 bits
                    chan    = bitand(bitshift(T3Record,-25),63);   % the next 6 bits
                    tcspc   = bitand(chan,15);
                    special = bitand(bitshift(T3Record,-31),1);    % the last bit

                    ind  = (special==1 & chan==63);

                    special = special.*chan;
                    
                otherwise
                    error('Illegal RecordType!');
            end;

            tmp  = sync(ind==1);
            tmp(tmp==0) = 1;
            sync(ind) = tmp;
            sync = sync + WRAPAROUND*cumsum(ind.*sync);
            
            sync(ind)    = [];
            tcspc(ind)   = [];
            special(ind) = [];
            chan(ind)    = [];
            loc = num - find(ind==0,1,'last');
                        
        end
        
        fclose(fid);
        
    else
        if (nargin<1)||isempty(name)
            fprintf(1,'\n\n      You have to specify a valid file-name. Aborted.\n');
        else
            fprintf(1,'\n\n      Could not open <%s>. Aborted.\n', name);
        end;
    end
end
