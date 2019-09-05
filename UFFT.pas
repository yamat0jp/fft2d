unit UFFT;

interface

uses System.SysUtils, System.Math;

type
  TFFTData = array of Extended;

  TFFT2dData = array of array of Extended;

  // 高速フーリエ変換(入力信号,出力実部,出力虚部)
procedure fft2d(n: integer; InRe2d: TFFT2dData; out OutRe, OutIm: TFFT2dData);

procedure fft(InRe: TFFTData; var OutRe, OutIm: TFFTData);
procedure _fft(InRe: TFFTData; var OutRe, OutIm: TFFTData);

// *******************窓関数*****************************
// ハン窓関数
procedure WinHanning(var data: TFFTData);
// ハミング窓関数
procedure WinHamming(var data: TFFTData);
// ガウス窓関数
procedure WinGauss(var data: TFFTData; m: integer = 1);
// ブラックマンハリス窓関数
procedure WinBlackmanHarris(var data: TFFTData);
// ブラックマンナットール窓関数
procedure WinBlackmanNuttall(var data: TFFTData);
// フラップトップ窓関数
procedure WinFlapTop(var data: TFFTData);
// 半正弦波窓関数
procedure WinHalfSin(var data: TFFTData);

implementation

procedure fft2d(n: integer; InRe2d: TFFT2dData; out OutRe, OutIm: TFFT2dData);
var
  s, Fr, Fi, Fr2, Fi2: TFFTData;
  j: integer;
  i: integer;
begin
  SetLength(s, n);
  SetLength(OutRe, n, n);
  SetLength(OutIm, n, n);
  for j := 0 to n - 1 do
  begin
    for i := 0 to n - 1 do
      s[i] := InRe2d[i, j];
    winHanning(s);
    fft(s, Fr, Fi);
    for i := 0 to n - 1 do
    begin
      OutRe[i, j] := Fr[i];
      OutIm[i, j] := Fi[i];
    end;
  end;
  for i := 0 to n - 1 do
  begin
    for j := 0 to n - 1 do
      s[j] := OutRe[i, j];
    winHanning(s);
    fft(s, Fr, Fi);
    for j := 0 to n - 1 do
      s[j] := OutIm[i, j];
    winHanning(s);
    fft(s, Fr2, Fi2);
    for j := 0 to n - 1 do
    begin
      OutRe[i, j] := Fr[j] - Fi2[j];
      OutIm[i, j] := Fi[j] + Fr2[j];
    end;
  end;
  Finalize(s);
  Finalize(Fr);
  Finalize(Fi);
  Finalize(Fr2);
  Finalize(Fi2);
end;

procedure fft(InRe: TFFTData; var OutRe, OutIm: TFFTData);
var
  i: integer;
  InN: integer; // 入力データ数
  n: integer; // 補正後データ数
begin
  InN := Length(InRe);       InN:=128;
  // データ数が2の乗数に満たない場合は0のデータを追加する
  i := 1;
  while InN >= Power(2, i) do
    inc(i);
  n := trunc(IntPower(2, i));
  if InN < n then
  begin
    SetLength(InRe, n);
    for i := InN to n - 1 do
      InRe[i] := 0;
  end;
  // 高速フーリエ変換
  _fft(InRe, OutRe, OutIm);
end;

procedure _fft(InRe: TFFTData; var OutRe, OutIm: TFFTData);
var
  n: integer;
  i: integer;
  InIm: TFFTData; // 複素数の虚部
  ct1, ct2, ct3: integer;
  TmpRe, TmpIm: Extended;
  nfft: array [0 .. 3] of integer;
  fcos, fsin: TFFTData;
  tmp: Extended;
  noblk: integer;
  cntb: array [0 .. 1] of integer;
begin
  n := Length(InRe);
  SetLength(InIm, n);
  for i := 0 to n - 1 do
    InIm[i] := 0;

  ct2 := 1;
  for ct1 := 1 to Length(InIm) - 2 do
  begin
    TmpRe := 0;
    TmpIm := 0;
    if ct1 < ct2 then
    begin
      TmpRe := InRe[ct1 - 1];
      InRe[ct1 - 1] := InRe[ct2 - 1];
      InRe[ct2 - 1] := TmpRe;
      TmpIm := InIm[ct1 - 1];
      InIm[ct1 - 1] := InIm[ct2 - 1];
      InIm[ct2 - 1] := TmpIm;
    end;
    ct3 := Length(InIm) div 2;
    while ct3 < ct2 do
    begin
      ct2 := ct2 - ct3;
      ct3 := ct3 div 2;
    end;
    ct2 := ct2 + ct3;
  end;

  nfft[0] := floor(Log2(Length(InIm)) / Log2(2));
  SetLength(fcos, n);
  SetLength(fsin, n);
  fcos[0] := 1;
  fsin[0] := 0;

  for ct1 := 1 to nfft[0] do
  begin
    nfft[2] := floor(System.Math.Power(2, ct1));
    nfft[1] := n div nfft[2];
    nfft[3] := nfft[2] div 2;
    for ct2 := 1 to nfft[3] do
    begin
      tmp := -Pi / nfft[3] * ct2;
      fcos[ct2] := cos(tmp);
      fsin[ct2] := sin(tmp);
    end;
    for ct2 := 1 to nfft[1] do
    begin
      noblk := nfft[2] * (ct2 - 1);
      for ct3 := 0 to nfft[3] - 1 do
      begin
        cntb[0] := noblk + ct3;
        cntb[1] := cntb[0] + nfft[3];
        TmpRe := InRe[cntb[1]] * fcos[ct3] - InIm[cntb[1]] * fsin[ct3];
        TmpIm := InIm[cntb[1]] * fcos[ct3] + InRe[cntb[1]] * fsin[ct3];
        InRe[cntb[1]] := InRe[cntb[0]] - TmpRe;
        InIm[cntb[1]] := InIm[cntb[0]] - TmpIm;
        InRe[cntb[0]] := InRe[cntb[0]] + TmpRe;
        InIm[cntb[0]] := InIm[cntb[0]] + TmpIm;
      end;
    end;
  end;
  SetLength(OutRe, Length(InRe) div 2);
  SetLength(OutIm, Length(InIm) div 2);
  for i := 0 to (n div 2) - 1 do
  begin
    OutRe[i] := InRe[i]; // 実部
    OutIm[i] := InIm[i]; // 虚部
  end;
end;

// ハン窓関数
procedure WinHanning(var data: TFFTData);
var
  i, n: integer;
begin
  n := Length(data);
  for i := 0 to n - 1 do
  begin
    data[i] := (0.5 - 0.5 * cos(2 * Pi * i / (n - 1))) * data[i];
  end;
end;

// ハミング窓関数
procedure WinHamming(var data: TFFTData);
var
  i, n: integer;
begin
  n := Length(data);
  for i := 0 to n - 1 do
  begin
    data[i] := (0.54 - 0.46 * cos(2 * Pi * i / (n - 1))) * data[i];
  end;
end;

// ガウス窓関数
procedure WinGauss(var data: TFFTData; m: integer = 1);
var
  i, n: integer;
begin
  n := Length(data);
  for i := 0 to n - 1 do
  begin
    data[i] := Exp(-2 * Power(m, 2) / Power(n - 1, 2) * Power(i - (n - 1) / 2,
      2)) * data[i];
  end;
end;

// ブラックマンハリス窓関数
procedure WinBlackmanHarris(var data: TFFTData);
var
  i, n: integer;
begin
  n := Length(data);
  for i := 0 to n - 1 do
  begin
    data[i] := (0.35875 - 0.48829 * cos(2 * Pi * i / (n - 1)) + 0.14128 *
      cos(4 * Pi * i / (n - 1)) - 0.01168 * cos(6 * Pi * i / (n - 1)))
      * data[i];
  end;
end;

// ブラックマンナットール窓関数
procedure WinBlackmanNuttall(var data: TFFTData);
var
  i, n: integer;
begin
  n := Length(data);
  for i := 0 to n - 1 do
  begin
    data[i] := (0.3635819 - 0.4891775 * cos(2 * Pi * i / (n - 1)) + 0.1365995 *
      cos(4 * Pi * i / (n - 1)) - 0.0106411 * cos(6 * Pi * i / (n - 1)))
      * data[i];
  end;
end;

// フラップトップ窓関数
procedure WinFlapTop(var data: TFFTData);
var
  i, n: integer;
begin
  n := Length(data);
  for i := 0 to n - 1 do
  begin
    data[i] := (1 - 1.930 * cos(2 * Pi * i / (n - 1)) + 1.290 *
      cos(4 * Pi * i / (n - 1)) - 0.388 * cos(6 * Pi * i / (n - 1)) + 0.032 *
      cos(8 * Pi * i / (n - 1))) * data[i];
  end;
end;

// 半正弦波窓関数
procedure WinHalfSin(var data: TFFTData);
var
  i, n: integer;
begin
  n := Length(data);
  for i := 0 to n - 1 do
  begin
    data[i] := sin(Pi * i / (n - 1)) * data[i];
  end;
end;

end.
