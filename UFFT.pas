unit UFFT;

interface
uses System.SysUtils,System.Math;
type
  TFFTData=array of Extended;

//高速フーリエ変換(入力信号,出力実部,出力虚部)
procedure fft(InRe:TFFTData;var OutRe,OutIm:TFFTData);
procedure _fft(InRe:TFFTData;var OutRe,OutIm:TFFTData);

//*******************窓関数*****************************
//ハン窓関数
procedure WinHanning(var data:TFFTData);
//ハミング窓関数
procedure WinHamming(var data:TFFTData);
//ガウス窓関数
procedure WinGauss(var data:TFFTData;m:integer=1);
//ブラックマンハリス窓関数
procedure WinBlackmanHarris(var data:TFFTData);
//ブラックマンナットール窓関数
procedure WinBlackmanNuttall(var data:TFFTData);
//フラップトップ窓関数
procedure WinFlapTop(var data:TFFTData);
//半正弦波窓関数
procedure WinHalfSin(var data:TFFTData);

implementation


procedure fft(InRe:TFFTData;var OutRe,OutIm:TFFTData);
var i  :Integer;
    InN:Integer;//入力データ数
    n  :Integer;//補正後データ数
begin
  InN:=Length(InRe);
  //データ数が2の乗数に満たない場合は0のデータを追加する
  i:=1;
  while InN > Power(2,i) do inc(i);
  n:=trunc(IntPower(2,i));
  if InN < n then
  begin
    setlength(InRe,n);
    for i := InN to n-1 do InRe[i]:=0;
  end;
  //高速フーリエ変換
  _fft(InRe,OutRe,OutIm);
end;

procedure _fft(InRe:TFFTData;var OutRe,OutIm:TFFTData);
var n:Integer;
    i:integer;
    InIm:TFFTData; //複素数の虚部
    ct1,ct2,ct3:integer;
    TmpRe,TmpIm:extended;
    nfft:array[0..3] of integer;
    fcos,fsin:TFFTData;
    tmp:extended;
    noblk:integer;
    cntb:array[0..1] of integer;
begin
  n:=Length(InRe);
  setlength(InIm,n);
  for i := 0 to n-1 do InIm[i]:=0;

  ct2:=1;
  for ct1 := 1 to length(InIm)-2 do
  begin
    TmpRe:=0;
    TmpIm:=0;
    if ct1<ct2 then
    begin
      TmpRe:=InRe[ct1-1];
      InRe[ct1-1]:=InRe[ct2-1];
      InRe[ct2-1]:=TmpRe;
      TmpIm:=InIm[ct1-1];
      InIm[ct1-1]:=InIm[ct2-1];
      InIm[ct2-1]:=TmpIm;
    end;
    ct3:=length(InIm) div 2;
    while ct3<ct2 do
    begin
      ct2:=ct2-ct3;
      ct3:=ct3 div 2;
    end;
    ct2:=ct2+ct3;
  end;

  nfft[0]:=floor(Log2(length(InIm))/Log2(2));
  SetLength(fcos,n);
  SetLength(fsin,n);
  fcos[0]:=1;
  fsin[0]:=0;

  for ct1 := 1 to nfft[0] do
  begin
    nfft[2]:=floor(System.math.Power(2,ct1));
    nfft[1]:=n div nfft[2];
    nfft[3]:=nfft[2] div 2;
    for ct2 := 1 to nfft[3] do
    begin
      tmp:=-Pi/nfft[3]*ct2;
      fcos[ct2]:=cos(tmp);
      fsin[ct2]:=sin(tmp);
    end;
    for ct2 := 1 to nfft[1] do
    begin
      noblk:=nfft[2]*(ct2-1);
      for ct3 := 0 to nfft[3]-1 do
      begin
        cntb[0]:=noblk+ct3;
        cntb[1]:=cntb[0]+nfft[3];
        TmpRe:=InRe[cntb[1]]*fcos[ct3]-InIm[cntb[1]]*fsin[ct3];
        TmpIm:=InIm[cntb[1]]*fcos[ct3]+InRe[cntb[1]]*fsin[ct3];
        InRe[cntb[1]]:=InRe[cntb[0]]-TmpRe;
        InIm[cntb[1]]:=InIm[cntb[0]]-TmpIm;
        InRe[cntb[0]]:=InRe[cntb[0]]+TmpRe;
        InIm[cntb[0]]:=InIm[cntb[0]]+TmpIm;
      end;
    end;
  end;
  setlength(OutRe,length(InRe) div 2);
  setlength(OutIm,length(InIm) div 2);
  for i := 0 to (n div 2)-1 do
  begin
    OutRe[i]:=InRe[i];  //実部
    OutIm[i]:=InIm[i];  //虚部
  end;
end;

//ハン窓関数
procedure WinHanning(var data:TFFTData);
var i,n:integer;
begin
  n:=length(data);
  for i := 0 to n-1 do
  begin
    data[i]:=(
               0.5 - 0.5*Cos(2*Pi*i/(n-1))
             )*data[i];
  end;
end;
//ハミング窓関数
procedure WinHamming(var data:TFFTData);
var i,n:integer;
begin
  n:=length(data);
  for i := 0 to n-1 do
  begin
    data[i]:=(
               0.54 - 0.46 * Cos(2*Pi*i/(n-1))
             )*data[i];
  end;
end;
//ガウス窓関数
procedure WinGauss(var data:TFFTData;m:integer=1);
var i,n:integer;
begin
  n:=length(data);
  for i := 0 to n-1 do
  begin
    data[i]:=Exp(
               -2 * power(m,2) / power(n-1,2) * power(i-(n-1)/2,2)
             )*data[i];
  end;
end;
//ブラックマンハリス窓関数
procedure WinBlackmanHarris(var data:TFFTData);
var i,n:integer;
begin
  n:=length(data);
  for i := 0 to n-1 do
  begin
    data[i]:=(0.35875-0.48829*cos(2*Pi*i/(n-1))
                     +0.14128*cos(4*Pi*i/(n-1))
                     -0.01168*cos(6*Pi*i/(n-1))
             )*data[i];
  end;
end;
//ブラックマンナットール窓関数
procedure WinBlackmanNuttall(var data:TFFTData);
var i,n:integer;
begin
  n:=length(data);
  for i := 0 to n-1 do
  begin
    data[i]:=(0.3635819-0.4891775*cos(2*Pi*i/(n-1))
                       +0.1365995*cos(4*Pi*i/(n-1))
                       -0.0106411*cos(6*Pi*i/(n-1))
             )*data[i];
  end;
end;
//フラップトップ窓関数
procedure WinFlapTop(var data:TFFTData);
var i,n:integer;
begin
  n:=length(data);
  for i := 0 to n-1 do
  begin
    data[i]:=(1-1.930*Cos(2*Pi*i/(n-1))
               +1.290*Cos(4*Pi*i/(n-1))
               -0.388*Cos(6*Pi*i/(n-1))
               +0.032*Cos(8*Pi*i/(n-1))
             )*data[i];
  end;
end;
//半正弦波窓関数
procedure WinHalfSin(var data:TFFTData);
var i,n:integer;
begin
  n:=length(data);
  for i := 0 to n-1 do
  begin
    data[i]:=Sin(pi*i/(n-1))*data[i];
  end;
end;

end.
