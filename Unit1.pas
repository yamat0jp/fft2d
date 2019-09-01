unit Unit1;

interface

uses
  System.SysUtils, System.Types, System.UITypes, System.Classes,
  System.Variants,
  FMX.Types, FMX.Controls, FMX.Forms, FMX.Dialogs, FMX.Ani, FMX.Layouts,
  FMX.Gestures,
  FMX.StdCtrls, FMX.ListBox, FMX.Controls.Presentation, FMX.Edit,
  FMX.TabControl,
  FMX.Objects;

type
  TForm1 = class(TForm)
    ToolbarHolder: TLayout;
    ToolbarPopup: TPopup;
    ToolbarPopupAnimation: TFloatAnimation;
    ToolBar1: TToolBar;
    ToolbarApplyButton: TButton;
    ToolbarCloseButton: TButton;
    ToolbarAddButton: TButton;
    OpenDialog1: TOpenDialog;
    TabControl1: TTabControl;
    TabItem1: TTabItem;
    TabItem2: TTabItem;
    Button1: TButton;
    Button2: TButton;
    Edit1: TEdit;
    Label1: TLabel;
    Button3: TButton;
    ListBox1: TListBox;
    Button4: TButton;
    Button5: TButton;
    Edit2: TEdit;
    Button6: TButton;
    Image1: TImage;
    Image2: TImage;
    Image3: TImage;
    procedure Button1Click(Sender: TObject);
    procedure Button6Click(Sender: TObject);
    procedure Button2Click(Sender: TObject);
    procedure FormCreate(Sender: TObject);
    procedure FormDestroy(Sender: TObject);
  private
    nx, ny: Integer;
    sr, si: array of array of Double;
    aPr, aPt: array of Double;
    meanC1, meanC2: Double;
    procedure fft2D;
    procedure calcPolarSpectra;
    procedure drawPolarSpectra;
    procedure sortingSmall(var a: array of Double; var ID: array of Integer;
      n: Integer);
    { private êÈåæ }
  public
    { public êÈåæ }
  end;

var
  Form1: TForm1;

implementation

{$R *.fmx}

uses FMX.Graphics, System.Math, Unit2, UFFT;

procedure main(n: Integer; sr, si: PDouble); stdcall;
  external 'CppMathLibrary.dll' name 'fft2d';

procedure sortingLarge(a: Pointer; ID: Pointer; n: Integer); stdcall;
  external 'CppMathLibrary.dll';

const
  maxRadius = 128;
  maxTheta = 180;

procedure TForm1.Button1Click(Sender: TObject);
begin
  if OpenDialog1.Execute = true then
  begin
    Image1.Bitmap.LoadFromFile(OpenDialog1.FileName);
    Edit1.Text := ExtractFileName(OpenDialog1.FileName);
  end;
end;

procedure TForm1.Button2Click(Sender: TObject);
begin
  fft2D;
  calcPolarSpectra;
  drawPolarSpectra;
end;

procedure TForm1.Button6Click(Sender: TObject);
const
  numEntry = 100;
  weightC = 1;
type
  TFeatureVector = record
    name: string;
    aPr: array [0 .. maxRadius - 1] of Double;
    aPt: array [0 .. maxTheta - 1] of Double;
    meanC1, meanC2: Double;
  end;
var
  i: Integer;
  k: Integer;
  distC, distR, distT: Double;
  a, b: Double;
  str: string;
  dist: array [0 .. numEntry - 1] of Double;
  ID: array [0 .. numEntry - 1] of Integer;
  fv: array [0 .. numEntry - 1] of TFeatureVector;
begin
  fft2D;
  calcPolarSpectra;
  for k := 0 to numEntry - 1 do
  begin
    distR := 0.0;
    distT := 0.0;
    distC := (fv[k].meanC1 - meanC1) * (fv[k].meanC1 - meanC1) +
      (fv[k].meanC2 - meanC2) * (fv[k].meanC2 - meanC2);
    distC := Sqrt(distC) * weightC;
    for i := 0 to maxRadius - 1 do
    begin
      a := fv[k].aPr[i];
      b := aPr[i];
      distR := distR + (a - b) * (a - b);
    end;
    distR := Sqrt(distR);
    for i := 0 to maxTheta - 1 do
    begin
      a := fv[k].aPt[i];
      b := aPt[i];
      distT := distT + (a - b) * (a - b);
    end;
    distT := Sqrt(distT);
    dist[k] := distR + distT + distC;
  end;
  for i := 0 to numEntry - 1 do
    ID[i] := i;
  sortingSmall(dist, ID, numEntry);
  for i := 1 to 5 do
  begin
    str := i.ToString + ' ' + ID[i].ToString + ' ' + fv[ID[i]].name + ' ' +
      dist[i].ToString;
    ListBox1.Items.Add(str)
  end;
end;

procedure TForm1.calcPolarSpectra;
var
  j: Integer;
  i: Integer;
  max: Double;
  rr, th: Integer;
  x, y: Double;
  Pr, Pt: array of Double;
begin
  SetLength(Pr, maxRadius);
  SetLength(Pt, maxTheta);
  for i := 0 to maxRadius - 1 do
    Pr[i] := 0.0;
  for i := 0 to maxTheta - 1 do
    Pt[i] := 0.0;
  for j := 0 to (ny div 2) + 10 do
    for i := 0 to (nx div 2) + 10 do
    begin
      x := i - nx / 2;
      y := ny / 2 - j;
      rr := Trunc(Sqrt(x * x + y * y) + 0.5);
      th := Trunc(ArcTan2(y, x) * 180.0 / pi + 0.5);
      if (rr <> 0) and (rr < maxRadius) and (th >= 0) and (th < maxTheta) then
      begin
        Pr[rr] := Pr[rr] + si[i, j];
        Pt[th] := Pt[th] + si[i, j];
      end;
    end;
  aPr[0] := (Pr[0] + Pr[1]) / 2.0;
  aPt[0] := (Pt[0] + Pt[1]) / 2.0;
  aPr[maxRadius - 1] := (Pr[maxRadius - 2] + Pr[maxRadius - 1]) / 2.0;
  aPt[maxTheta - 1] := (Pt[maxTheta - 2] + Pt[maxTheta - 1]) / 2.0;
  for i := 0 to maxRadius - 1 do
    aPr[i] := (Pr[i - 1] + Pr[i] + Pr[i + 1]) / 3.0;
  for i := 0 to maxTheta - 1 do
    aPt[i] := (Pt[i - 1] + Pt[i] + Pt[i + 1]) / 3.0;
  max := 0.0;
  for i := 0 to maxRadius - 1 do
    if aPr[i] > max then
      max := aPr[i];
  for i := 0 to maxRadius - 1 do
    aPr[i] := aPr[i] / max;
  max := 0.0;
  for i := 0 to maxTheta - 1 do
    if aPt[i] > max then
      max := aPt[i];
  for i := 0 to maxTheta - 1 do
    aPt[i] := aPt[i] / max;
  Finalize(sr);
  Finalize(si);
end;

procedure TForm1.drawPolarSpectra;
var
  i, h, x0, y0: Integer;
begin
  Image3.Bitmap.Assign(Image2.Bitmap);
  Image3.Bitmap.Width := 290;
  Image3.Bitmap.Height := 170;
  Image3.Bitmap.Clear(TAlphaColors.White);
  h := Image3.Bitmap.Height - 50;
  x0 := 20;
  y0 := Image3.Bitmap.Height - 20;
  with Image3.Bitmap.Canvas do
  begin
    BeginScene;
    Fill.Color := TAlphaColors.Blue;
    FillText(RectF(0, 0, 20, 20), 'p(r)', false, 1, [], TTextAlign.Center,
      TTextAlign.Center);
    DrawLine(PointF(x0, y0), PointF(x0, maxRadius + 10), 1);
    DrawLine(PointF(x0, y0 + 5), PointF(x0, y0 - h - 10), 1);
    DrawLine(PointF(x0 + maxRadius, y0), PointF(x0 + maxRadius, y0 + 5), 1);
    DrawLine(PointF(x0 + maxRadius / 2, y0), PointF(x0 + maxRadius / 2,
      y0 + 5), 1);
    for i := 1 to maxRadius - 1 do
      DrawLine(PointF(x0 + i - 1, y0 - aPr[i - 1] * h),
        PointF(x0 + i, y0 - aPr[i] * h), 1);
    h := Image3.Bitmap.Height - 50;
    x0 := 60 + maxRadius;
    y0 := Image3.Bitmap.Height - 20;
    FillText(RectF(x0, 0, x0 + 20, 20), 'q(É∆)', false, 1, [],
      TTextAlign.Center);
    DrawLine(PointF(x0, y0), PointF(x0 + maxTheta / 2 + 10, y0), 1);
    DrawLine(PointF(x0, y0 - 5), PointF(x0, y0 - h - 10), 1);
    DrawLine(PointF(x0, y0), PointF(x0, y0 - aPt[0] * h), 1);
    DrawLine(PointF(x0 + maxTheta / 2, y0), PointF(x0 + maxTheta, y0 + 5), 1);
    DrawLine(PointF(x0 + maxTheta / 4, y0), PointF(x0 + maxTheta / 4,
      y0 + 5), 1);
    for i := 1 to maxTheta div 2 - 1 do
      DrawLine(PointF(x0 + i - 1, y0 - aPt[2 * (i - 1)] * h),
        PointF(x0 + i, y0 - aPt[2 * i] * h), 1);
    EndScene;
  end;
end;

procedure TForm1.fft2D;
const
  GMAX = 255;
  dRange = 60;
var
  j: Integer;
  i: Integer;
  gray: Integer;
  rr, gg, bb: Double;
  light: Double;
  db: Double;
  a, b, c, max: Double;
  data: TBitmapData;
  col: TAlphaColorRec;
  fft: Tfft;
begin
  nx := Image1.Bitmap.Width;
  ny := Image1.Bitmap.Height;
  SetLength(sr, nx, ny);
  SetLength(si, nx, ny);
  Image1.Bitmap.Map(TMapAccess.ReadWrite, data);
  try
    for j := 0 to ny - 1 do
      for i := 0 to nx - 1 do
      begin
        col.Color := data.GetPixel(i, j);
        rr := col.R;
        gg := col.G;
        bb := col.b;
        light := 0.299 * rr + 0.587 * gg + 0.114 * bb;
        meanC1 := meanC1 + rr - light;
        meanC2 := meanC2 + bb - light;
        gray := Trunc(light);
        col.R := gray;
        col.G := gray;
        col.b := gray;
        data.SetPixel(i, j, col.Color);
        sr[i, j] := light;
        si[i, j] := 0.0;
      end;
  finally
    Image1.Bitmap.Unmap(data);
  end;
  meanC1 := meanC1 / (nx * ny);
  meanC2 := meanC2 / (nx * ny);
  fft := Tfft.Create;
  try
    fft.fft2D(nx, PDouble(sr), PDouble(si));
  finally
    fft.Free;
  end;
  max := 0.0;
  for j := 0 to ny - 1 do
    for i := 0 to nx - 1 do
    begin
      a := sr[i, j] / GMAX;
      b := si[i, j] / GMAX;
      c := a * a + b * b;
      sr[i, j] := c;
      if c = 0.0 then
        db := 0.0
      else
        db := 10.0 * logN(10, c) + dRange;
      if db < 0.0 then
        db := 0.0;
      if max < db then
        max := db;
      si[i, j] := db;
    end;
  for j := 0 to (ny div 2) - 1 do
  begin
    for i := 0 to (nx div 2) - 1 do
    begin
      a := sr[i, j];
      sr[i, j] := sr[i + nx div 2, j + ny div 2];
      sr[i + nx div 2, j + ny div 2] := a;
      a := si[i, j];
      si[i, j] := si[i + nx div 2, j + ny div 2];
      si[i + nx div 2, j + ny div 2] := a;
    end;
    for i := nx div 2 to nx - 1 do
    begin
      a := sr[i, j];
      sr[i, j] := sr[i - nx div 2, j + ny div 2];
      sr[i - nx div 2, j + ny div 2] := a;
      a := si[i, j];
      si[i, j] := si[i - nx div 2, j + ny div 2];
      si[i - nx div 2, j + ny div 2] := a;
    end;
  end;
  Image2.Bitmap.Assign(Image1.Bitmap);
  Image2.Bitmap.Clear(TAlphaColors.Black);
  Image2.Bitmap.Map(TMapAccess.Write, data);
  try
    for j := 0 to ny - 1 do
      for i := 0 to nx - 1 do
      begin
        gray := Trunc(255.0 * si[i, j] / max);
        col.R := gray;
        col.G := gray;
        col.b := gray;
        data.SetPixel(i, j, col.Color);
      end;
  finally
    Image2.Bitmap.Unmap(data);
  end;
end;

procedure TForm1.FormCreate(Sender: TObject);
begin
  SetLength(aPr, maxRadius);
  SetLength(aPt, maxTheta);

end;

procedure TForm1.FormDestroy(Sender: TObject);
begin
  Finalize(aPr);
  Finalize(aPt);
end;

procedure TForm1.sortingSmall(var a: array of Double; var ID: array of Integer;
  n: Integer);
var
  k: Integer;
  k1: Integer;
  min: Double;
  sid: Integer;
begin
  for k := 0 to n - 1 do
  begin
    min := a[k];
    sid := ID[k];
    for k1 := k + 1 to n do
      if min > a[k1] then
      begin
        a[k] := a[k1];
        a[k1] := min;
        min := a[k];
        ID[k] := ID[k1];
        ID[k1] := sid;
        sid := ID[k];
      end;
  end;
end;

end.
