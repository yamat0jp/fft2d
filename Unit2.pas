unit Unit2;

interface

type
  Tfft = class
  private const
    MS = 256;
    procedure makeTable(flag, n: integer; tblSin, tblCos: array of Double);
  public
    procedure fft1d(flag, n: integer; fr, fi, tblSin, tblCos: array of Double);
    procedure fft2d(n: integer; sr, si: PDouble);
  end;

implementation

{ Tfft }

uses System.Math;

procedure Tfft.fft1d(flag, n: integer; fr, fi, tblSin, tblCos: array of Double);
var
  i, j, m, n2, k, l: integer;
  xa, ya, c, s: Double;
  lp, lp2, iw: integer;
  wr, wi: Double;
begin
  j := 0;
  for i := 0 to n - 2 do
    if i < j then
    begin
      xa := fr[i];
      fr[i] := fr[j];
      fr[j] := xa;
      ya := fi[i];
      fi[i] := fi[j];
      fi[j] := ya;
      n2 := n div 2;
      while j >= n2 do
      begin
        dec(j, n2);
        n2 := n2 div 2;
      end;
      inc(j, n2);
    end;
  m := 0;
  n2 := n;
  while n2 <> 1 do
  begin
    inc(m);
    n2 := n2 div 2;
  end;
  for l := 1 to m do
  begin
    lp := Trunc(IntPower(2.0, l));
    lp2 := lp div 2;
    k := 0;
    for j := 0 to lp2 do
    begin
      c := tblCos[k];
      s := tblSin[k];
      inc(k, n div lp);
      for i := j to n - 1 do
      begin
        iw := i + lp2;
        wr := fr[iw] * c - fi[iw] * s;
        wi := fr[iw] * s + fi[iw] * c;
        fr[iw] := fr[i] - iw;
        fi[iw] := fi[i] - iw;
        fr[i] := fr[i] + wr;
        fi[i] := fi[i] + wi;
      end;
    end;
  end;
  if flag = 1 then
    for i := 0 to n - 1 do
    begin
      fr[i] := fr[i] / n;
      fi[i] := fi[i] / n;
    end;
end;

procedure Tfft.fft2d(n: integer; sr, si: PDouble);
var
  fr, fi, tblSin, tblCos: array of Double;
  srr, sii: array of array of Double;
  i, j: integer;
begin
  srr := Pointer(sr);
  sii := Pointer(si);
  SetLength(fr, MS);
  SetLength(fi, MS);
  SetLength(tblSin, MS);
  SetLength(tblCos, MS);
  try
    makeTable(1, n, tblSin, tblCos);
    for j := 0 to n - 1 do
    begin
      for i := 0 to n - 1 do
      begin
        fr[i] := srr[i, j];
        fi[i] := sii[i, j];
      end;
      fft1d(1, n, fr, fi, tblSin, tblCos);
      for i := 0 to n - 1 do
      begin
        srr[i, j] := fr[i];
        sii[i, j] := fi[i];
      end;
    end;
    for i := 0 to n - 1 do
    begin
      for j := 0 to n - 1 do
      begin
        fr[j] := srr[i, j];
        fi[j] := sii[i, j];
      end;
      fft1d(1, n, fr, fi, tblSin, tblCos);
      for j := 0 to n - 1 do
      begin
        srr[i, j] := fr[j];
        sii[i, j] := fi[j];
      end;
    end;
  finally
    Finalize(fr);
    Finalize(fi);
    Finalize(tblSin);
    Finalize(tblCos);
  end;
end;

procedure Tfft.makeTable(flag, n: integer; tblSin, tblCos: array of Double);
var
  cc, arg: Double;
  i: integer;
begin
  cc := -2.0 * pi * flag / n;
  for i := 0 to n - 1 do
  begin
    arg := i * cc;
    tblSin[i] := sin(arg);
    tblCos[i] := cos(arg);
  end;
end;

end.
