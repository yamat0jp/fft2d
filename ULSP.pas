unit ULSP;

interface
uses System.SysUtils,System.Math;

type
  TLSPDataSingleton=record
    PastSec:Extended;//経過時間(秒)
    Data:Extended;   //入力信号
  end;
  TLSPData=array of TLSPDataSingleton;
  TLSPSpecSingleton=record
    Hz:Extended;
    Data:Extended;
  end;
  TLSPSpec=array of TLSPSpecSingleton;

procedure lsp(InData:TLSPData;var OutData:TLSPSpec;
              SHz,EHz,PHz:Extended);

implementation

uses Unit1;

//データ,出力データ,解析する最小周波数,解析する最大周波数,解析するステップ周波数
procedure lsp(InData:TLSPData;var OutData:TLSPSpec;
              SHz,EHz,PHz:Extended);
var i,j:integer;
    imax:integer;
    c,s,ys1,yc1,cc1,cc2,ss1,ss2,cs1,cs2:Extended;
    d,p,p_u,p_d:Extended;
    w,wt:extended;
    DataCt:integer;
    DataSum,DataAve:extended;
    sigma:Extended;
    coswt,sinwt:Extended;
begin
  DataCt:=Length(InData);

  DataSum:=0;
  for i := 0 to DataCt-1 do
  begin
    DataSum:=DataSum+InData[i].Data;
  end;
  DataAve:=DataSum/DataCt;

  sigma:=0;
  for i := 0 to DataCt-1 do
  begin
    sigma:=sigma+power(InData[i].Data-DataAve,2);
  end;

  imax:=floor((EHz-SHz) / PHz);
  for i := 0 to imax do
  begin
    c:=0;
    s:=0;
    ys1:=0;
    yc1:=0;
    cc2:=0;
    ss2:=0;
    cs2:=0;
    w:=2*pi*(SHz+PHz*i);
    for j := 0 to DataCt-1 do
    begin
      wt:=w*InData[j].PastSec;
      coswt:=cos(wt);
      sinwt:=sin(wt);
      c:=c+coswt;
      s:=s+sinwt;
      yc1:=yc1+(InData[j].Data-DataAve)*coswt;
      ys1:=ys1+(InData[j].Data-DataAve)*sinwt;
      cc2:=cc2+power(coswt,2);
      ss2:=ss2+power(sinwt,2);
      cs2:=cs2+coswt*sinwt;
    end;
    cc1:=cc2-c*c;
    ss1:=ss2-s*s;
    cs1:=cs2-c*s;
    d:=cc1*ss1-cs1*cs1;
    p_u:=ss1*yc1*yc1+cc1*ys1*ys1-2*cs1*yc1*ys1;
    p_d:=sigma*d;
    if p_d=0 then p_d:=1;
    p:=p_u/p_d;
    setlength(OutData,i+1);
    OutData[i].Hz:=SHz+PHz*i;
    OutData[i].Data:=sqrt(abs(p));
  end;
end;

end.
