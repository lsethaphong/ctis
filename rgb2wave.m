function Maplambda = rgb2wave(Im);
% 5 AUG 04 by L. Sethaphong
%                    proportion of total intensity
% index lamda	R	G	B  
% 1     400	131	0	181
% 2     410	126	0	219
% 3     420	106	0	255
% 4     430	61	0	255
% 5     440	0	0	255
% 6     450	0	70	255
% 7     460	0	123	255
% 8     470	0	169	255
% 9     480	0	213	255
% 10    490	0	255	255
% 11    500	0	255	135
% 12    510	0	255	0
% 13    520	54	255	0
% 14    530	94	255	0
% 15    540	129	255	0
% 16    550	163	255	0
% 17    560	195	255	0
% 18    570	225	255	0
% 19    580	255	255	0
% 20    590	255	223	0
% 21    600	255	190	0
% 22    610	255	155	0
% 23    620	255	119	0
% 24    630	255	79	0
% 25    640	255	33	0
% 26    650	255	0	0
% 27    660	255	0	0
% 28    670	255	0	0
% 29    680	255	0	0
% 30    690	255	0	0
% 31    700	255	0	0
% 32    710   230 0   0
% 33    720   210 0   0
% 34    730   190 0   0
% 35    740   170 0   0
% 36    750   150 0   0
% 37    760   130 0   0
% 38    770   110 0   0
% 39    780   90  0   0
% 40    790   70  0   0
% 41    800   50  0   0
% voxel displacement
% create storage array
%try
RES = size(Im);
mapIm = zeros(RES(1),RES(2));
for i=1:RES(1) % rows
    for j=1:RES(2) % columns
        R = double(Im(i,j,1));
        G = double(Im(i,j,2));
        B = double(Im(i,j,3));
        % floor lowest value
        % contribution of luminosity
     if (R+G+B) > 0
        pR = R/(R+G+B);
        pG = G/(R+G+B);
        pB = B/(R+G+B);
%
% 24    630	255	79	0
% 25    640	255	33	0
% 26    650	255	0	0 --
% 27    660	251	0	0
% 28    670	247	0	0
% 29    680	243	0	0
% 30    690	239	0	0
% 31    700	234	0	0
% 32    710   230 0   0
% 33    720   210 0   0
% 34    730   190 0   0
% 35    740   170 0   0
% 36    750   150 0   0
% 37    760   130 0   0
% 38    770   110 0   0
% 39    780   90  0   0
% 40    790   70  0   0
% 41    800   50  0   0        
% 630nm-800nm
        if pR > 0.75
            if pG > 0.20
               mapIm(i,j)=24;
            elseif pG > 0.10
               mapIm(i,j)=25;
            elseif R > 254
               mapIm(i,j)=26;
            elseif R > 250
               mapIm(i,j)=27;
            elseif R > 246
               mapIm(i,j)=28;
            elseif R > 242
               mapIm(i,j)=29;
            elseif R > 238
               mapIm(i,j)=30;
            elseif R > 234
               mapIm(i,j)=31;
            elseif R > 229
               mapIm(i,j)=32;                
            elseif R > 209
               mapIm(i,j)=33;                
            elseif R > 189
               mapIm(i,j)=34;                
            elseif R > 169
               mapIm(i,j)=35;                
            elseif R > 149
               mapIm(i,j)=36;                
            elseif R > 129
               mapIm(i,j)=37;                
            elseif R > 109
               mapIm(i,j)=38;                
            elseif R > 89
               mapIm(i,j)=39;                
            elseif R > 59
               mapIm(i,j)=40;                
            else
               mapIm(i,j)=41;
            end
%
% 12    510	0	255	0
% 13    520	54	255	0
% 14    530	94	255	0
% 15    540	129	255	0
% 16    550	163	255	0
% 17    560	195	255	0
% 18    570	225	255	0
% 19    580	255	255	0
% 20    590	255	223	0
% 21    600	255	190	0
% 22    610	255	155	0
% 23    620	255	119	0
% 510nm-630nm
        elseif (pR+pG) > 0.75
            ppG = pG/(pR+pG);
            if ppG > 0.9
                mapIm(i,j)=12;
            elseif ppG > 0.8 
                mapIm(i,j)=13;
            elseif ppG > 0.7
                mapIm(i,j)=14;
            elseif ppG > 0.65
                mapIm(i,j)=15;
            elseif ppG > 0.6
                mapIm(i,j)=16;
            elseif ppG > 0.55
                mapIm(i,j)=17;
            elseif ppG > 0.5
                mapIm(i,j)=18;
            elseif ppG > 0.47
                mapIm(i,j)=19;
            elseif ppG > 0.43
                mapIm(i,j)=20;
            elseif ppG > 0.40
                mapIm(i,j)=21;
            elseif ppG > 0.35
                mapIm(i,j)=22;
            else
                mapIm(i,j)=23;
            end
%
% 5     440	0	0	255
% 6     450	0	70	255
% 7     460	0	123	255
% 8     470	0	169	255
% 9     480	0	213	255
% 10    490	0	255	255
% 11    500	0	255	135        
% 440nm-500nm
        elseif (pG+pB) > 0.85
            ppB=pB/(pG+pB);
            if ppB > 0.95 
                mapIm(i,j)= 5;
            elseif ppB > 0.77
                mapIm(i,j)= 6;
            elseif ppB > 0.66
                mapIm(i,j)= 7;
            elseif ppB > 0.6 
                mapIm(i,j)= 8;
            elseif ppB > 0.5
                mapIm(i,j)= 9;
            elseif ppB > 0.4
                mapIm(i,j)= 10;
            else
                mapIm(i,j)= 11;
            end
%
% 1     400	131	0	181
% 2     410	126	0	219
% 3     420	106	0	255
% 4     430	61	0	255
% 400nm-430nm
        elseif (pR+pB) > 0.85
            ppB=pB/(pR+pB);
            if ppB > 0.71
                mapIm(i,j)= 4;
            elseif ppB > 0.69
                mapIm(i,j)= 3;
            elseif ppB > 0.59
                mapIm(i,j)= 2;
            else
                mapIm(i,j)= 1;
            end
        end
            
        % mapping
    end
    end
end
Maplambda = mapIm;
%catch
%    i
%    j
%MapLambda = zeros(RES,RES);
%end