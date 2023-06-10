using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Diagnostics;
using System.Diagnostics.CodeAnalysis;
using System.Globalization;
using System.Linq;
using System.Runtime.CompilerServices;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using System.Windows.Media.Media3D;

// ReSharper disable TailRecursiveCall
// ReSharper disable UnusedMember.Local

namespace Lab_1
{
    public partial class MainWindow
    {
        private const int Size = 600;
        private static int Offset => Size / 2;

        private double _degree = Math.PI / 4;

        private Point _startPoint = new Point(0, 0);

        private readonly Random _rnd = new Random();

        private static List<Point> _borderPoints;
        private static List<Point> _innerPoints;

        public event PropertyChangedEventHandler PropertyChanged;

        private void NotifyPropertyChanged([CallerMemberName] string propertyName = "")
        {
            PropertyChanged?.Invoke(this, new PropertyChangedEventArgs(propertyName));
        }

        public double Degree
        {
            get => _degree;
            set
            {
                if (Math.Abs(_degree - value) < 0.01)
                    return;

                _degree = value;
                RedrawGraph();
                TextBoxDegree.Content = Degree.ToString(CultureInfo.CurrentCulture);
                NotifyPropertyChanged();
            }
        }

        public MainWindow()
        {
            InitializeComponent();

            _borderPoints = GetRandomPolygonPoints(Offset - 40, 9);
            _innerPoints = GetRandomPolygonPoints(Offset - 130);
            RedrawGraph();
            Slider.TickFrequency = Math.PI / 64;
            Slider.Maximum = Math.PI * 2;
        }

        private void DrawAxis()
        {
            var geometry = new StreamGeometry
            {
                FillRule = FillRule.Nonzero
            };

            using (var ctx = geometry.Open())
            {
                ctx.BeginFigure(new Point(0, Offset), true, false);
                ctx.LineTo(new Point(Size, Offset), true, false);
                ctx.LineTo(new Point(Offset, 0), false, false);
                ctx.LineTo(new Point(Offset, Size), true, false);
            }

            AddGeometry(geometry);
        }

        private void CalcCoordinates(ref double x, ref double y)
        {
            var x0 = _startPoint.X;
            var y0 = _startPoint.Y;
            var cosDegree = Math.Cos(Degree);
            var sinDegree = Math.Sin(Degree);
            var fakeX = x0 + (x - x0) * cosDegree - (y - y0) * sinDegree;
            y = y0 + (x - x0) * sinDegree + (y - y0) * cosDegree;
            x = fakeX;
        }

        // ReSharper disable once UnusedMember.Local
        private void DrawMathFuncLine()
        {
            var geometry = new StreamGeometry
            {
                FillRule = FillRule.Nonzero
            };

            using (var ctx = geometry.Open())
            {
                double Y(double x) => x + 10;
                var init = true;
                for (var x = -Offset; x < Offset; x += 20)
                {
                    double fakeX = x;
                    var fakeY = Offset - Y(x + Offset);
                    CalcCoordinates(ref fakeX, ref fakeY);
                    if (fakeY < -Offset || fakeY > Offset || fakeX < -Offset || fakeX > Offset)
                        continue;

                    var point = new Point(Offset + fakeX, Offset + fakeY);
                    if (init)
                    {
                        ctx.BeginFigure(point, true, false);
                        init = false;
                    }

                    ctx.LineTo(point, true, false);
                }
            }

            AddGeometry(geometry);
        }

        private void RedrawGraph()
        {
            ClearGeometry();

            //Lab 1
            //DrawAxis();
            //DrawMathFuncLine();

            //Lab 2
            //var rotatedPolygonPoints = ReCalcPoints(_innerPoints);
            //var clippedPolygonPoints = SutherlandHodgman(rotatedPolygonPoints, _borderPoints);
            //DrawPolygon(clippedPolygonPoints, true);
            //DrawPolygon(_borderPoints, true);

            //Lab 3
            //Draw(500, 0, 120, Math.PI / 2);

            //Lab 4
            //_borderPoints = GetRandomPolygonPoints(Offset - 50, 100);
            //DrawPolygon(_borderPoints, true);

            //Lab 5
            //StreamGeometry.Clear();
            //DrawPolygon(_borderPoints);
            //SplitPolygon(_borderPoints);

            //Lab 6
            //DrawAxis();
            //DrawProjection();

            //Lab 7
            DrawCube(_vecDegree[1], _vecDegree[0]);

            //Lab 8
            //DrawLightPhong();
            //ImageForBitMap.Visibility = Visibility.Hidden;
        }

        private WriteableBitmap _writeableBitmap;

        private Vector _sourceLight = new Vector(300, -300, 300);

        private const int RadiusCircle = 125;
        private const int RadiusSquared = RadiusCircle * RadiusCircle;
        private const int SizeImage = 600;
        private const int OffsetImage = SizeImage / 2;

        [SuppressMessage("ReSharper.DPA", "DPA0002: Excessive memory allocations in SOH",
            MessageId = "type: Vector; size: 1155MB")]
        private void DrawLightPhong()
        {
            var first = false;
            if (_writeableBitmap == null)
            {
                _writeableBitmap = new WriteableBitmap(
                    SizeImage,
                    SizeImage,
                    96,
                    96,
                    PixelFormats.Bgr32,
                    null);
                RenderOptions.SetBitmapScalingMode(ImageForBitMap, BitmapScalingMode.NearestNeighbor);
                RenderOptions.SetEdgeMode(ImageForBitMap, EdgeMode.Aliased);
                ImageForBitMap.Stretch = Stretch.None;
                ImageForBitMap.HorizontalAlignment = HorizontalAlignment.Left;
                ImageForBitMap.VerticalAlignment = VerticalAlignment.Top;

                ImageForBitMap.Source = _writeableBitmap;
                first = true;
            }

            var viewVector = new Vector(0, 0, -250).Normalized(); // Направление обзора наблюдателя
            var lightDirection = _sourceLight.Normalized(); // Нормализованное направление источника света

            var ambientColor = new Vector(255, 255, 0); // Цвет фонового освещения
            var diffuseColor = new[] // Цвет диффузного освещения
            {
                new Vector(0, 255, 0),
                new Vector(0, 0, 255)
            };
            var specularColor = new Vector(255, 255, 0); // Цвет зеркального отражения
            const double materialShininess = 50; // Коэффициент зеркального отражения материала

            const double ambientIntensity = 0.2f; // Интенсивность фонового освещения
            const double diffuseIntensity = 0.8f; // Интенсивность диффузного освещения
            const double specularIntensity = 0.5f; // Интенсивность зеркального отражения

            var centerX = new[] {75, -75};
            var centerY = new[] {75, -75};

            var ambient = ambientIntensity * ambientColor;

            try
            {
                _writeableBitmap.Lock();
                if (first)
                {
                    for (var x = 0; x < SizeImage; x++)
                    {
                        for (var y = 0; y < SizeImage; y++)
                        {
                            DrawPixel(x, y, (int) ambientColor.X, (int) ambientColor.Y, (int) ambientColor.Z);
                        }
                    }
                }

                for (var i = 0; i < centerX.Length; i++)
                {
                    for (var x = -RadiusCircle + centerX[i]; x <= RadiusCircle + centerX[i]; x++)
                    {
                        var limit = Math.Sqrt(RadiusSquared - Math.Pow(x - centerX[i], 2));
                        for (var y = (int) -limit + centerY[i]; y <= limit + centerY[i]; y++)
                        {
                            var z = Math.Sqrt(RadiusSquared - Math.Pow(x - centerX[i], 2) -
                                              Math.Pow(y - centerY[i], 2));

                            var normalVector = new Vector((double) (x + centerX[i]) / RadiusCircle,
                                (double) (y + centerY[i]) / RadiusCircle, z / RadiusCircle).Normalized();

                            var diffuse = diffuseIntensity * diffuseColor[i] *
                                          Math.Max(Vector.Dot(normalVector, lightDirection), 0);
                            var reflection = Vector.Reflect(-lightDirection, normalVector).Normalized();
                            var specular = specularIntensity * specularColor *
                                           Math.Pow(Math.Max(Vector.Dot(reflection, viewVector), 0), materialShininess);

                            var light = ambient + diffuse + specular;
                            if (light.X > 254 || light.Y > 254 || light.Z > 254)
                                light = light.Normalized() * 255;

                            DrawPixel((x + centerX[i]) + OffsetImage, (y + centerY[i]) + OffsetImage, (int) light.X,
                                (int) light.Y, (int) light.Z);
                        }
                    }
                }
            }
            finally
            {
                _writeableBitmap.Unlock();
            }
        }

        private void DrawLightGyro()
        {
            var first = false;
            if (_writeableBitmap == null)
            {
                _writeableBitmap = new WriteableBitmap(
                    SizeImage,
                    SizeImage,
                    96,
                    96,
                    PixelFormats.Bgr32,
                    null);
                RenderOptions.SetBitmapScalingMode(ImageForBitMap, BitmapScalingMode.NearestNeighbor);
                RenderOptions.SetEdgeMode(ImageForBitMap, EdgeMode.Aliased);
                ImageForBitMap.Stretch = Stretch.None;
                ImageForBitMap.HorizontalAlignment = HorizontalAlignment.Left;
                ImageForBitMap.VerticalAlignment = VerticalAlignment.Top;

                ImageForBitMap.Source = _writeableBitmap;
                first = true;
            }

            new Vector(0, 0, -250).Normalized();
            var lightDirection = _sourceLight.Normalized(); // Нормализованное направление источника света

            var ambientColor = new Vector(255, 255, 0); // Цвет фонового освещения
            var diffuseColor = new[] // Цвет диффузного освещения
            {
                new Vector(0, 255, 0),
                new Vector(0, 0, 255)
            };

            const double diffuseIntensity = 0.8f; // Интенсивность диффузного освещения

            var centerX = new[] {75, -75};
            var centerY = new[] {75, -75};

            try
            {
                _writeableBitmap.Lock();
                if (first)
                {
                    for (var x = 0; x < SizeImage; x++)
                    {
                        for (var y = 0; y < SizeImage; y++)
                        {
                            DrawPixel(x, y, (int) ambientColor.X, (int) ambientColor.Y, (int) ambientColor.Z);
                        }
                    }
                }

                for (var i = 0; i < centerX.Length; i++)
                {
                    for (var x = -RadiusCircle + centerX[i]; x <= RadiusCircle + centerX[i]; x++)
                    {
                        var limit = Math.Sqrt(RadiusSquared - Math.Pow(x - centerX[i], 2));
                        for (var y = (int) -limit + centerY[i]; y <= limit + centerY[i]; y++)
                        {
                            var z = Math.Sqrt(RadiusSquared - Math.Pow(x - centerX[i], 2) -
                                              Math.Pow(y - centerY[i], 2));

                            var normalVector = new Vector((double) (x + centerX[i]) / RadiusCircle,
                                (double) (y + centerY[i]) / RadiusCircle, z / RadiusCircle).Normalized();

                            var diffuse = diffuseIntensity * diffuseColor[i] *
                                          Math.Max(Vector.Dot(normalVector, lightDirection), 0);

                            if (diffuse.X > 254 || diffuse.Y > 254 || diffuse.Z > 254)
                                diffuse = diffuse.Normalized() * 255;

                            DrawPixel((x + centerX[i]) + OffsetImage, (y + centerY[i]) + OffsetImage, (int) diffuse.X,
                                (int) diffuse.Y, (int) diffuse.Z);
                        }
                    }
                }
            }
            finally
            {
                _writeableBitmap.Unlock();
            }
        }

        private class Vector
        {
            public readonly double X;
            public readonly double Y;
            public readonly double Z;

            public Vector(double x, double y, double z)
            {
                X = x;
                Y = y;
                Z = z;
            }

            private double Length()
                => Math.Sqrt(X * X + Y * Y + Z * Z);

            public Vector Normalized()
            {
                var length = Length();
                return length > 0 ? new Vector(X / length, Y / length, Z / length) : new Vector(0, 0, 0);
            }

            public static double Dot(Vector v1, Vector v2)
                => v1.X * v2.X + v1.Y * v2.Y + v1.Z * v2.Z;

            public static Vector Reflect(Vector v, Vector normal)
            {
                var dotProduct = Dot(v, normal);
                return new Vector(v.X - 2 * dotProduct * normal.X, 
                    v.Y - 2 * dotProduct * normal.Y, 
                    v.Z - 2 * dotProduct * normal.Z);
            }

            public static Vector operator *(Vector vector, double k)
                => new Vector(vector.X * k, vector.Y * k, vector.Z * k);

            public static Vector operator *(double k, Vector vector)
                => vector * k;

            public static Vector operator +(Vector v1, Vector v2)
                => new Vector(v1.X + v2.X, v1.Y + v2.Y, v1.Z + v2.Z);

            public static Vector operator -(Vector v)
                => new Vector(-v.X, -v.Y, -v.Z);

            public static Vector operator -(Vector vectorA, Vector vectorB) 
                => new Vector(vectorA.X - vectorB.X, vectorA.Y - vectorB.Y, vectorA.Z - vectorB.Z);

            public static Vector Cross(Vector vectorA, Vector vectorB) =>
                new Vector(vectorA.Y * vectorB.Z - vectorA.Z * vectorB.Y,
                    vectorA.Z * vectorB.X - vectorA.X * vectorB.Z, 
                    vectorA.X * vectorB.Y - vectorA.Y * vectorB.X);
        }

        //writeableBitmap should already be locked
        private unsafe void DrawPixel(int x, int y, int r, int g, int b)
        {
            var pBackBuffer = _writeableBitmap.BackBuffer;
            pBackBuffer += x * _writeableBitmap.BackBufferStride;
            pBackBuffer += y * 4;
            var colorData = r << 16; // R
            colorData |= g << 8; // G
            colorData |= b << 0; // B
            *(int*) pBackBuffer = colorData;
            _writeableBitmap.AddDirtyRect(new Int32Rect(x, y, 1, 1));
        }

        private void ImageForBitMap_OnMouseMove(object sender, MouseEventArgs e)
        {
            ProcessMouseInput(e);
        }

        private readonly double[] _vecDegree = {0, 0};

        private void ProcessMouseInput(MouseEventArgs e)
        {
            if (e.RightButton != MouseButtonState.Pressed && e.LeftButton != MouseButtonState.Pressed)
                return;

            var point = e.GetPosition(ImageForBitMap);
            double x = (int) point.Y - OffsetImage;
            double y = (int) point.X - OffsetImage;
            if (e.LeftButton == MouseButtonState.Pressed)
            {
                const double z = 125;
                _sourceLight = new Vector(x / OffsetImage, y / OffsetImage, z / OffsetImage);
            }

            if (e.RightButton == MouseButtonState.Pressed)
            {
                _vecDegree[0] = x / OffsetImage * Math.PI;
                _vecDegree[1] = y / OffsetImage * Math.PI;
            }

            RedrawGraph();
        }

        private void DrawProjection()
        {
            const double z1 = 700;
            const double z2 = 500;
            const double d = z1 / z2;
            var points = new List<Point>
            {
                new Point(100, 100),
                new Point(100, -100),
                new Point(-100, -100),
                new Point(-100, 100),
                new Point(100, 100)
            };
            var n = points.Count;
            for (var i = 0; i < n; i++)
            {
                points.Add(new Point(d * points[i].X, d * points[i].Y));
            }

            DrawPolygon(points);

            var geometry = new StreamGeometry
            {
                FillRule = FillRule.Nonzero
            };
            using (var ctx = geometry.Open())
            {
                ctx.BeginFigure(points[0] /*decoy*/, false, false);
                for (var i = 0; i < n - 1; i++)
                {
                    ctx.LineTo(new Point(points[i].X + Offset, points[i].Y + Offset), false, false);
                    ctx.LineTo(new Point(points[i + 5].X + Offset, points[i + 5].Y + Offset), true, false);
                }
            }

            AddGeometry(geometry);
        }

        private void SplitPolygon(List<Point> borderPoints)
        {
            var t = (int) (Degree / 0.025);
            var list = new List<Point>(borderPoints);
            //StreamGeometry.Clear();
            using (var ctx = StreamGeometry.Open())
            {
                ctx.BeginFigure(new Point(0, 0), false, false);
                var n = list.Count;
                var i = 0;
                while (n > 3)
                {
                    i = (i + 1) % n;
                    var lPoint = i != 0 ? list[i - 1] : list[n - 1];
                    var rPoint = i != n - 1 ? list[i + 1] : list[0];
                    var point = list[i];
                    if (CheckOrientation(lPoint, rPoint, point) > 0)
                        continue;

                    list.Remove(list[i]);
                    if (n >= t)
                    {
                        ctx.LineTo(new Point(lPoint.X + Offset, lPoint.Y + Offset), false, false);
                        ctx.LineTo(new Point(rPoint.X + Offset, rPoint.Y + Offset), true, false);
                    }

                    i += 1;
                    n = list.Count;
                }

                if (t != 0)
                    return;
                
                StreamGeometry.Clear();
                ctx.BeginFigure(new Point(list[0].X + Offset, list[0].Y + Offset), true, true);
                ctx.LineTo(new Point(list[1].X + Offset, list[1].Y + Offset), true, true);
                ctx.LineTo(new Point(list[2].X + Offset, list[2].Y + Offset), true, true);
            }
        }

        private static double CheckOrientation(Point a, Point b, Point c) // >0 - C слева. <0 - C справа. =0 - C коллинеарна
        {
            return (b.X - a.X) * (c.Y - a.Y) - (c.X - a.X) * (b.Y - a.Y);
        }

        private void Draw(double x1, double y1, double lineLength, double angle)
        {
            var x2 = x1 + lineLength * Math.Cos(angle);
            var y2 = y1 + lineLength * Math.Sin(angle);

            var angle1 = angle + Math.PI / 2;
            var x3 = x2 + lineLength * Math.Cos(angle1);
            var y3 = y2 + lineLength * Math.Sin(angle1);

            var angle2 = angle + Math.PI / 2;
            var x4 = x1 + lineLength * Math.Cos(angle2);
            var y4 = y1 + lineLength * Math.Sin(angle2);

            var fi = Degree;
            var angle3 = angle - fi;
            var x5 = x3 + lineLength * Math.Cos(fi) * Math.Cos(angle3);
            var y5 = y3 + lineLength * Math.Cos(fi) * Math.Sin(angle3);

            var geometry = new StreamGeometry
            {
                FillRule = FillRule.Nonzero
            };

            using (var ctx = geometry.Open())
            {
                ctx.BeginFigure(new Point((int) Math.Round(x1), (int) Math.Round(Size - y1)), false, false);
                ctx.LineTo(new Point((int) Math.Round(x2), (int) Math.Round(Size - y2)), true, false);
                ctx.LineTo(new Point((int) Math.Round(x3), (int) Math.Round(Size - y3)), true, false);
                ctx.LineTo(new Point((int) Math.Round(x4), (int) Math.Round(Size - y4)), true, false);
                ctx.LineTo(new Point((int) Math.Round(x1), (int) Math.Round(Size - y1)), true, false);
            }

            AddGeometry(geometry);

            lineLength *= 0.71;
            if (lineLength <= 4)
                return;

            Draw(x2, y2, lineLength, angle - fi);
            Draw(x5, y5, lineLength, angle + fi);
        }

        private void AddGeometry(Geometry g)
        {
            GeometryGroup.Children.Add(g);
        }

        private void ClearGeometry()
        {
            while (GeometryGroup.Children.Count > 0)
                GeometryGroup.Children.Remove(GeometryGroup.Children[0]);
        }

        private void TextBoxBase_OnTextChanged(object sender, TextChangedEventArgs e)
        {
            var split = ((TextBox) sender).Text.Split(' ');
            if (split.Length < 2)
                return;

            if (!double.TryParse(split[0], out var x))
                return;

            if (!double.TryParse(split[1], out var y))
                return;

            _startPoint = new Point(x, y);
            RedrawGraph();
        }

        private List<Point> GetRandomPolygonPoints(double radius, int n = 0)
        {
            var points = new List<Point>();
            if (n <= 2)
            {
                for (var a = 0.0; a < 2.0 * Math.PI;) // full circle
                {
                    var x = Offset + radius * Math.Cos(a);
                    var y = Offset + radius * Math.Sin(a);
                    points.Add(new Point(x, y));

                    a += (20.0 + 80.0 * _rnd.NextDouble()) * Math.PI / 180.0;
                }
            }
            else
            {
                var sumN = Math.PI / n;
                for (var a = 0.0; a < 2.0 * Math.PI; a += sumN) // full circle
                {
                    var bufRadius = radius * (_rnd.NextDouble() % 0.9 + 0.1);
                    var x = Offset + bufRadius * Math.Cos(a);
                    var y = Offset + bufRadius * Math.Sin(a);
                    points.Add(new Point(x, y));
                }
            }

            return points;
        }

        private List<Point> ReCalcPoints(List<Point> points)
        {
            points = new List<Point>(points);
            for (var i = 0; i < points.Count; i++)
            {
                var x = points[i].X;
                var y = points[i].Y;
                CalcCoordinates(ref x, ref y);
                points[i] = new Point(x, y);
            }

            return points;
        }

        private void DrawPolygon(List<Point> points, bool isShifted = false)
        {
            if (points.Count == 0)
                return;

            var geometry = new StreamGeometry
            {
                FillRule = FillRule.Nonzero
            };

            using (var ctx = geometry.Open())
            {
                if (isShifted)
                {
                    ctx.BeginFigure(points[0], true, true);
                    for (var i = 1; i < points.Count; i++)
                        ctx.LineTo(points[i], true, false);
                }
                else
                {
                    ctx.BeginFigure(new Point(points[0].X + Offset, points[0].Y + Offset), true, true);
                    for (var i = 1; i < points.Count; i++)
                        ctx.LineTo(new Point(points[i].X + Offset, points[i].Y + Offset), true, false);
                }
            }

            geometry.Freeze();
            AddGeometry(geometry);
        }

        public static List<Point> SutherlandHodgman(List<Point> subjectPolygon, List<Point> clipPolygon)
        {
            var resultPolygon = subjectPolygon;

            foreach (var clipEdgeStart in clipPolygon)
            {
                var clipEdgeEndIndex = clipPolygon.IndexOf(clipEdgeStart) == clipPolygon.Count - 1
                    ? 0
                    : clipPolygon.IndexOf(clipEdgeStart) + 1;
                var clipEdgeEnd = clipPolygon[clipEdgeEndIndex];

                var inputList = resultPolygon;
                resultPolygon = new List<Point>();
                if (inputList.Count == 0)
                    break;

                var s = inputList[inputList.Count - 1];

                foreach (var e in inputList)
                {
                    if (IsInside(e, clipEdgeStart, clipEdgeEnd))
                    {
                        if (!IsInside(s, clipEdgeStart, clipEdgeEnd))
                        {
                            var point = ComputeIntersection(s, e, clipEdgeStart, clipEdgeEnd);
                            if (point != null)
                                resultPolygon.Add(point.Value);
                        }

                        resultPolygon.Add(e);
                    }
                    else if (IsInside(s, clipEdgeStart, clipEdgeEnd))
                    {
                        var point = ComputeIntersection(s, e, clipEdgeStart, clipEdgeEnd);
                        if (point != null)
                            resultPolygon.Add(point.Value);
                    }

                    s = e;
                }
            }

            return resultPolygon;
        }

        private static bool IsInside(Point testPoint, Point clipEdgeStart, Point clipEdgeEnd)
        {
            return (clipEdgeEnd.X - clipEdgeStart.X) * (testPoint.Y - clipEdgeStart.Y) >
                   (clipEdgeEnd.Y - clipEdgeStart.Y) * (testPoint.X - clipEdgeStart.X);
        }

        private static Point? ComputeIntersection(Point line1Start, Point line1End, Point line2Start, Point line2End)
        {
            var denominator = (line2End.Y - line2Start.Y) * (line1End.X - line1Start.X) -
                              (line2End.X - line2Start.X) * (line1End.Y - line1Start.Y);

            if (denominator == 0) // Lines are parallel
                return null;

            var ua = ((line2End.X - line2Start.X) * (line1Start.Y - line2Start.Y) -
                      (line2End.Y - line2Start.Y) * (line1Start.X - line2Start.X)) / denominator;

            var intersectionX = line1Start.X + ua * (line1End.X - line1Start.X);
            var intersectionY = line1Start.Y + ua * (line1End.Y - line1Start.Y);

            return new Point(intersectionX, intersectionY);
        }
        
        public class Plane
        {
            public double A { get; }
            public double B { get; }
            public double C { get; }
            public double D { get; }

            public Plane(double x1, double y1, double z1, double x2, double y2, double z2, double x3, double y3, double z3)
            {
                A = (y2 - y1) * (z3 - z1) - (z2 - z1) * (y3 - y1);
                B = (z2 - z1) * (x3 - x1) - (x2 - x1) * (z3 - z1);
                C = (x2 - x1) * (y3 - y1) - (y2 - y1) * (x3 - x1);
                D = -(A * x1 + B * y1 + C * z1);
            }

            public double CalculateZ(double x, double y)
            {
                return (-A * x - B * y - D) / C;
            }
        }
        
        private void DrawCube(double degreeX, double degreeY)
        {
            if (_writeableBitmap == null)
            {
                _writeableBitmap = new WriteableBitmap(
                    SizeImage,
                    SizeImage,
                    96,
                    96,
                    PixelFormats.Bgr32,
                    null);
                RenderOptions.SetBitmapScalingMode(ImageForBitMap, BitmapScalingMode.NearestNeighbor);
                RenderOptions.SetEdgeMode(ImageForBitMap, EdgeMode.Aliased);
                ImageForBitMap.Stretch = Stretch.None;
                ImageForBitMap.HorizontalAlignment = HorizontalAlignment.Left;
                ImageForBitMap.VerticalAlignment = VerticalAlignment.Top;

                ImageForBitMap.Source = _writeableBitmap;
            }
            
            var viewVector = new Vector(0, 0, -250).Normalized();
            var lightDirection = _sourceLight.Normalized();

            var ambientColor = new Vector(200, 200, 75);
            var diffuseColor = new Vector(0, 255, 0);
            var specularColor = new Vector(255, 255, 0);
            const double materialShininess = 50;

            const double ambientIntensity = 0.2f;
            const double diffuseIntensity = 0.8f;
            const double specularIntensity = 0.5f;

            var ambient = ambientIntensity * ambientColor;

            const int defSize = 200;
            var size = (int) Math.Round(defSize * mashtab);
            var halfSize = size / 2;
            var centerX = _writeableBitmap.PixelWidth / 2;
            var centerY = _writeableBitmap.PixelHeight / 2;

            _writeableBitmap.Lock();
            for (var x = 0; x < _writeableBitmap.PixelWidth; x++)
            for (var y = 0; y < _writeableBitmap.PixelHeight; y++)
                DrawPixel(x, y, (int) ambientColor.X, (int) ambientColor.Y, (int) ambientColor.Z);

            var vertices = new[]
            {
                new Point3D(-halfSize, -halfSize, -halfSize),
                new Point3D(halfSize, -halfSize, -halfSize),
                new Point3D(halfSize, halfSize, -halfSize),
                new Point3D(-halfSize, halfSize, -halfSize),
                new Point3D(-halfSize, -halfSize, halfSize),
                new Point3D(halfSize, -halfSize, halfSize),
                new Point3D(halfSize, halfSize, halfSize),
                new Point3D(-halfSize, halfSize, halfSize)
            };

            var depthBuffer = new double[_writeableBitmap.PixelWidth, _writeableBitmap.PixelHeight];

            var faces = new[]
            {
                new[] {0, 1, 2, 3},
                new[] {1, 5, 6, 2},
                new[] {4, 5, 6, 7},
                new[] {0, 4, 7, 3},
                new[] {0, 1, 5, 4},
                new[] {3, 2, 6, 7}
            };

            foreach (var face in faces)
            {
                var screenVertices = new Vector[face.Length];
                for (var i = 0; i < face.Length; i++)
                {
                    var vertex = vertices[face[i]];

                    var rotatedY = vertex.Y * Math.Cos(degreeX) - vertex.Z * Math.Sin(degreeX);
                    var rotatedZ = vertex.Y * Math.Sin(degreeX) + vertex.Z * Math.Cos(degreeX);

                    var rotatedX = vertex.X * Math.Cos(degreeY) + rotatedZ * Math.Sin(degreeY);
                    
                    var screenX = centerX + rotatedX;
                    var screenY = centerY + rotatedY;

                    screenVertices[i] = new Vector(screenX, screenY, rotatedZ);
                }

                var plane = new Plane(screenVertices[0].X, screenVertices[0].Y, screenVertices[0].Z,
                    screenVertices[1].X, screenVertices[1].Y, screenVertices[1].Z,
                    screenVertices[2].X, screenVertices[2].Y, screenVertices[2].Z);

                var normal = new Vector(
                    (screenVertices[0].X + screenVertices[1].X + screenVertices[2].X + screenVertices[3].X) / 4,
                    (screenVertices[0].Y + screenVertices[1].Y + screenVertices[2].Y + screenVertices[3].Y) / 4,
                    (screenVertices[0].Z + screenVertices[1].Z + screenVertices[2].Z + screenVertices[3].Z) / 4)
                    .Normalized();

                for (var y = 0; y < _writeableBitmap.PixelHeight; y++)
                {
                    var scanlineMinX = int.MaxValue;
                    var scanlineMaxX = int.MinValue;

                    for (var i = 0; i < screenVertices.Length; i++)
                    {
                        var currentVertex = screenVertices[i];
                        var nextVertex = screenVertices[(i + 1) % screenVertices.Length];

                        if ((!(currentVertex.Y <= y) || !(nextVertex.Y > y)) && (!(nextVertex.Y <= y) || !(currentVertex.Y > y)))
                            continue;
                        
                        var xIntersection =
                            (nextVertex.X - currentVertex.X) * (y - currentVertex.Y) /
                            (nextVertex.Y - currentVertex.Y) + currentVertex.X;

                        if (xIntersection < scanlineMinX)
                            scanlineMinX = (int) xIntersection;
                        
                        if (xIntersection > scanlineMaxX)
                            scanlineMaxX = (int) xIntersection;
                    }

                    if (scanlineMinX > scanlineMaxX)
                        continue;

                    for (var x = scanlineMinX; x <= scanlineMaxX; x++)
                    {
                        var depth = plane.CalculateZ(x,y);
                        if (!(depth > depthBuffer[x, y]) && depthBuffer[x, y] != 0)
                            continue;
                        
                        var diffuse = diffuseIntensity * diffuseColor *
                                      Math.Max(Vector.Dot(normal, lightDirection), 0);
                        var reflection = Vector.Reflect(-lightDirection, normal).Normalized();
                        var specular = specularIntensity * specularColor *
                                       Math.Pow(Math.Max(Vector.Dot(reflection, viewVector), 0), materialShininess);

                        var light = ambient + diffuse + specular;
                        if (light.X > 254 || light.Y > 254 || light.Z > 254)
                            light = light.Normalized() * 255;
                        DrawPixel(x, y, (int) light.X, (int) light.Y, (int) light.Z);
                        
                        depthBuffer[x, y] = depth;
                    }
                }
            }

            _writeableBitmap.Unlock();
        }

        private void DrawLine(int x1, int y1, int x2, int y2, int r, int g, int b)
        {
            var dx = Math.Abs(x2 - x1);
            var dy = Math.Abs(y2 - y1);
            var sx = x1 < x2 ? 1 : -1;
            var sy = y1 < y2 ? 1 : -1;
            var err = dx - dy;

            while (true)
            {
                if (x1 >= 0 && x1 < _writeableBitmap.PixelWidth && y1 >= 0 && y1 < _writeableBitmap.PixelHeight)
                {
                    DrawPixel(x1, y1, r, g, b);
                }

                if (x1 == x2 && y1 == y2)
                {
                    break;
                }

                var e2 = 2 * err;

                if (e2 > -dy)
                {
                    err -= dy;
                    x1 += sx;
                }

                if (e2 < dx)
                {
                    err += dx;
                    y1 += sy;
                }
            }
        }

        private double mashtab = 1.0;
        private void ImageForBitMap_OnMouseWheel(object sender, MouseWheelEventArgs e)
        {
            var delta = e.Delta > 0 ? 0.01 : -0.01;
            mashtab += delta;
            if (mashtab < 0.3)
                mashtab = 0.3;
            else if (mashtab > 1.5)
                mashtab = 1.5;
            RedrawGraph();
        }
    }
}