+++
date = "2017-12-20T21:04:13-08:00"
title = "Pygments"
description = "A set of code samples, to demonstrate code formatting by Pygments."
categories = [ "Demo" ]
tags = [ "Pygments", "reStructuredText", "Hugo", "test" ]

+++

Pygments
########

.. class:: sidebar narrow

.. contents::

Code examples from
`Rosetta Code <http://rosettacode.org/wiki/Mandelbrot_set>`__
used to provide samples of code formatting by
`Pygments <http://pygments.org/>`__.


Awk
===

.. code:: Awk

   BEGIN {
     XSize=59; YSize=21;
     MinIm=-1.0; MaxIm=1.0;MinRe=-2.0; MaxRe=1.0;
     StepX=(MaxRe-MinRe)/XSize; StepY=(MaxIm-MinIm)/YSize;
     for(y=0;y<YSize;y++)
     {
       Im=MinIm+StepY*y;
       for(x=0;x<XSize;x++)
           {
         Re=MinRe+StepX*x; Zr=Re; Zi=Im;
         for(n=0;n<30;n++)
             {
           a=Zr*Zr; b=Zi*Zi;
           if(a+b>4.0) break;
           Zi=2*Zr*Zi+Im; Zr=a-b+Re;
         }
         printf "%c",62-n;
       }
       print "";
     }
     exit;
   }


C
=

.. code:: C

   /* 
   c program:
   --------------------------------
    1. draws Mandelbrot set for Fc(z)=z*z +c
    using Mandelbrot algorithm ( boolean escape time )
   -------------------------------         
   2. technique of creating ppm file is  based on the code of Claudio Rocchini
   http://en.wikipedia.org/wiki/Image:Color_complex_plot.jpg
   create 24 bit color graphic file ,  portable pixmap file = PPM 
   see http://en.wikipedia.org/wiki/Portable_pixmap
   to see the file use external application ( graphic viewer)
    */
   #include <stdio.h>
   #include <math.h>
   int main()
   {
            /* screen ( integer) coordinate */
          int iX,iY;
          const int iXmax = 800; 
          const int iYmax = 800;
          /* world ( double) coordinate = parameter plane*/
          double Cx,Cy;
          const double CxMin=-2.5;
          const double CxMax=1.5;
          const double CyMin=-2.0;
          const double CyMax=2.0;
          /* */
          double PixelWidth=(CxMax-CxMin)/iXmax;
          double PixelHeight=(CyMax-CyMin)/iYmax;
          /* color component ( R or G or B) is coded from 0 to 255 */
          /* it is 24 bit color RGB file */
          const int MaxColorComponentValue=255; 
          FILE * fp;
          char *filename="new1.ppm";
          char *comment="# ";/* comment should start with # */
          static unsigned char color[3];
          /* Z=Zx+Zy*i  ;   Z0 = 0 */
          double Zx, Zy;
          double Zx2, Zy2; /* Zx2=Zx*Zx;  Zy2=Zy*Zy  */
          /*  */
          int Iteration;
          const int IterationMax=200;
          /* bail-out value , radius of circle ;  */
          const double EscapeRadius=2;
          double ER2=EscapeRadius*EscapeRadius;
          /*create new file,give it a name and open it in binary mode  */
          fp= fopen(filename,"wb"); /* b -  binary mode */
          /*write ASCII header to the file*/
          fprintf(fp,"P6\n %s\n %d\n %d\n %d\n",comment,iXmax,iYmax,MaxColorComponentValue);
          /* compute and write image data bytes to the file*/
          for(iY=0;iY<iYmax;iY++)
          {
               Cy=CyMin + iY*PixelHeight;
               if (fabs(Cy)< PixelHeight/2) Cy=0.0; /* Main antenna */
               for(iX=0;iX<iXmax;iX++)
               {         
                          Cx=CxMin + iX*PixelWidth;
                          /* initial value of orbit = critical point Z= 0 */
                          Zx=0.0;
                          Zy=0.0;
                          Zx2=Zx*Zx;
                          Zy2=Zy*Zy;
                          /* */
                          for (Iteration=0;Iteration<IterationMax && ((Zx2+Zy2)<ER2);Iteration++)
                          {
                              Zy=2*Zx*Zy + Cy;
                              Zx=Zx2-Zy2 +Cx;
                              Zx2=Zx*Zx;
                              Zy2=Zy*Zy;
                          };
                          /* compute  pixel color (24 bit = 3 bytes) */
                          if (Iteration==IterationMax)
                          { /*  interior of Mandelbrot set = black */
                             color[0]=0;
                             color[1]=0;
                             color[2]=0;                           
                          }
                       else 
                          { /* exterior of Mandelbrot set = white */
                               color[0]=255; /* Red*/
                               color[1]=255;  /* Green */ 
                               color[2]=255;/* Blue */
                          };
                          /*write color to the file*/
                          fwrite(color,1,3,fp);
                  }
          }
          fclose(fp);
          return 0;
   }


Obfuscated C
============

.. code:: C

   main(k){float i,j,r,x,y=-16;while(puts(""),y++<15)for(x
   =0;x++<84;putchar(" .:-;!/>)|&IH%*#"[k&15]))for(i=k=r=0;
   j=r*r-i*i-2+x/25,i=2*r*i+y/10,j*j+i*i<11&&k++<111;r=j);}


C++
===

.. code:: C++

   #include <cstdlib>
   #include <complex>
    
   // get dimensions for arrays
   template<typename ElementType, std::size_t dim1, std::size_t dim2>
    std::size_t get_first_dimension(ElementType (&a)[dim1][dim2])
   {
     return dim1;
   }
    
   template<typename ElementType, std::size_t dim1, std::size_t dim2>
    std::size_t get_second_dimension(ElementType (&a)[dim1][dim2])
   {
     return dim2;
   }
    
    
   template<typename ColorType, typename ImageType>
    void draw_Mandelbrot(ImageType& image,                                   //where to draw the image
                         ColorType set_color, ColorType non_set_color,       //which colors to use for set/non-set points
                         double cxmin, double cxmax, double cymin, double cymax,//the rect to draw in the complex plane
                         unsigned int max_iterations)                          //the maximum number of iterations
   {
     std::size_t const ixsize = get_first_dimension(image);
     std::size_t const iysize = get_first_dimension(image);
     for (std::size_t ix = 0; ix < ixsize; ++ix)
       for (std::size_t iy = 0; iy < iysize; ++iy)
       {
         std::complex<double> c(cxmin + ix/(ixsize-1.0)*(cxmax-cxmin), cymin + iy/(iysize-1.0)*(cymax-cymin));
         std::complex<double> z = 0;
         unsigned int iterations;
    
         for (iterations = 0; iterations < max_iterations && std::abs(z) < 2.0; ++iterations) 
           z = z*z + c;
    
         image[ix][iy] = (iterations == max_iterations) ? set_color : non_set_color;
    
       }
   }


Clojure
=======

.. code:: Clojure

   (ns mandelbrot
     (:refer-clojure :exclude [+ * <])
     (:use (clojure.contrib complex-numbers)
           (clojure.contrib.generic [arithmetic :only [+ *]]
                                    [comparison :only [<]]
                                    [math-functions :only [abs]])))
   (defn mandelbrot? [z]
     (loop [c 1
            m (iterate #(+ z (* % %)) 0)]
       (if (and (> 20 c)
                (< (abs (first m)) 2) )
         (recur (inc c)
                (rest m))
         (if (= 20 c) true false))))
    
   (defn mandelbrot []
     (for [y (range 1 -1 -0.05)
   	x (range -2 0.5 0.0315)] 
       (if (mandelbrot? (complex x y)) "#" " ")))
    
   (println (interpose \newline (map #(apply str %) (partition 80 (mandelbrot)))))
    

Common Lisp
===========

.. code:: Common-Lisp

   (defpackage #:mandelbrot
     (:use #:cl))
    
   (in-package #:mandelbrot)
    
   (deftype pixel () '(unsigned-byte 8))
   (deftype image () '(array pixel))
    
   (defun write-pgm (image filespec)
     (declare (image image))
     (with-open-file (s filespec :direction :output :element-type 'pixel :if-exists :supersede)
       (let* ((width  (array-dimension image 1))
              (height (array-dimension image 0))
              (header (format nil "P5~A~D ~D~A255~A" #\Newline width height #\Newline #\Newline)))
         (loop for c across header
               do (write-byte (char-code c) s))
         (dotimes (row height)
           (dotimes (col width)
             (write-byte (aref image row col) s))))))
    
   (defparameter *x-max* 800)
   (defparameter *y-max* 800)
   (defparameter *cx-min* -2.5)
   (defparameter *cx-max* 1.5)
   (defparameter *cy-min* -2.0)
   (defparameter *cy-max* 2.0)
   (defparameter *escape-radius* 2)
   (defparameter *iteration-max* 40)
    
   (defun mandelbrot (filespec)
     (let ((pixel-width  (/ (- *cx-max* *cx-min*) *x-max*))
           (pixel-height (/ (- *cy-max* *cy-min*) *y-max*))
           (image (make-array (list *y-max* *x-max*) :element-type 'pixel :initial-element 0)))
       (loop for y from 0 below *y-max*
             for cy from *cy-min* by pixel-height
             do (loop for x from 0 below *x-max*
                      for cx from *cx-min* by pixel-width
                      for iteration = (loop with c = (complex cx cy)
                                            for iteration from 0 below *iteration-max*
                                            for z = c then (+ (* z z) c)
                                            while (< (abs z) *escape-radius*)
                                            finally (return iteration))
                      for pixel = (round (* 255 (/ (- *iteration-max* iteration) *iteration-max*)))
                      do (setf (aref image y x) pixel)))
       (write-pgm image filespec)))

Erlang
======

.. code:: Erlang

   -module(mandelbrot).
    
   -export([test/0]).
    
   magnitude(Z) ->
     R = complex:real(Z),
     I = complex:imaginary(Z),
     R * R + I * I.
    
   mandelbrot(A, MaxI, Z, I) ->
       case (I < MaxI) and (magnitude(Z) < 2.0) of
           true ->
               NZ = complex:add(complex:mult(Z, Z), A),
               mandelbrot(A, MaxI, NZ, I + 1);
           false ->
               case I of 
                   MaxI ->
                       $*;
                   _ ->
                       $ 
               end
       end.
    
   test() ->
       lists:map(
           fun(S) -> io:format("~s",[S]) end, 
           [
               [
                   begin 
                       Z = complex:make(X, Y),
                       mandelbrot(Z, 50, Z, 1)
                   end
               || X <- seq_float(-2, 0.5, 0.0315)
               ] ++ "\n"
           || Y <- seq_float(-1,1, 0.05)
           ] ),
       ok.
    
   % **************************************************
   % Copied from https://gist.github.com/andruby/241489
   % **************************************************
    
   seq_float(Min, Max, Inc, Counter, Acc) when (Counter*Inc + Min) >= Max -> 
     lists:reverse([Max|Acc]);
   seq_float(Min, Max, Inc, Counter, Acc) -> 
     seq_float(Min, Max, Inc, Counter+1, [Inc * Counter + Min|Acc]).
   seq_float(Min, Max, Inc) -> 
     seq_float(Min, Max, Inc, 0, []).
    
   % **************************************************
    

FORTRAN
=======

.. code:: FORTRAN

   program mandelbrot
    
     implicit none
     integer  , parameter :: rk       = selected_real_kind (9, 99)
     integer  , parameter :: i_max    =  800
     integer  , parameter :: j_max    =  600
     integer  , parameter :: n_max    =  100
     real (rk), parameter :: x_centre = -0.5_rk
     real (rk), parameter :: y_centre =  0.0_rk
     real (rk), parameter :: width    =  4.0_rk
     real (rk), parameter :: height   =  3.0_rk
     real (rk), parameter :: dx_di    =   width / i_max
     real (rk), parameter :: dy_dj    = -height / j_max
     real (rk), parameter :: x_offset = x_centre - 0.5_rk * (i_max + 1) * dx_di
     real (rk), parameter :: y_offset = y_centre - 0.5_rk * (j_max + 1) * dy_dj
     integer, dimension (i_max, j_max) :: image
     integer   :: i
     integer   :: j
     integer   :: n
     real (rk) :: x
     real (rk) :: y
     real (rk) :: x_0
     real (rk) :: y_0
     real (rk) :: x_sqr
     real (rk) :: y_sqr
    
     do j = 1, j_max
       y_0 = y_offset + dy_dj * j
       do i = 1, i_max
         x_0 = x_offset + dx_di * i
         x = 0.0_rk
         y = 0.0_rk
         n = 0
         do
           x_sqr = x ** 2
           y_sqr = y ** 2
           if (x_sqr + y_sqr > 4.0_rk) then
             image (i, j) = 255
             exit
           end if
           if (n == n_max) then
             image (i, j) = 0
             exit
           end if
           y = y_0 + 2.0_rk * x * y
           x = x_0 + x_sqr - y_sqr
           n = n + 1
         end do
       end do
     end do
     open  (10, file = 'out.pgm')
     write (10, '(a/ i0, 1x, i0/ i0)') 'P2', i_max, j_max, 255
     write (10, '(i0)') image
     close (10)
    
   end program mandelbrot


Go
==

.. code:: Go

   package main
    
   import "fmt"
   import "math/cmplx"
    
   func mandelbrot(a complex128) (z complex128) {
       for i := 0; i < 50; i++ {
           z = z*z + a
       }
       return
   }
    
   func main() {
       for y := 1.0; y >= -1.0; y -= 0.05 {
           for x := -2.0; x <= 0.5; x += 0.0315 {
               if cmplx.Abs(mandelbrot(complex(x, y))) < 2 {
                   fmt.Print("*")
               } else {
                   fmt.Print(" ")
               }
           }
           fmt.Println("")
       }
   }
   


Haskell
=======

.. code:: Haskell

   import Data.Complex
    
   mandelbrot a = iterate (\z -> z^2 + a) 0 !! 50
    
   main = mapM_ putStrLn [[if magnitude (mandelbrot (x :+ y)) < 2 then '*' else ' '
                              | x <- [-2, -1.9685 .. 0.5]]
                          | y <- [1, 0.95 .. -1]]
   


Java
====

.. code:: Java

   import java.awt.Graphics;
   import java.awt.image.BufferedImage;
   import javax.swing.JFrame;
    
   public class Mandelbrot extends JFrame {
    
       private final int MAX_ITER = 570;
       private final double ZOOM = 150;
       private BufferedImage I;
       private double zx, zy, cX, cY, tmp;
    
       public Mandelbrot() {
           super("Mandelbrot Set");
           setBounds(100, 100, 800, 600);
           setResizable(false);
           setDefaultCloseOperation(EXIT_ON_CLOSE);
           I = new BufferedImage(getWidth(), getHeight(), BufferedImage.TYPE_INT_RGB);
           for (int y = 0; y < getHeight(); y++) {
               for (int x = 0; x < getWidth(); x++) {
                   zx = zy = 0;
                   cX = (x - 400) / ZOOM;
                   cY = (y - 300) / ZOOM;
                   int iter = MAX_ITER;
                   while (zx * zx + zy * zy < 4 && iter > 0) {
                       tmp = zx * zx - zy * zy + cX;
                       zy = 2.0 * zx * zy + cY;
                       zx = tmp;
                       iter--;
                   }
                   I.setRGB(x, y, iter | (iter << 8));
               }
           }
       }
    
       @Override
       public void paint(Graphics g) {
           g.drawImage(I, 0, 0, this);
       }
    
       public static void main(String[] args) {
           new Mandelbrot().setVisible(true);
       }
   }


JavaScript
==========

.. code:: JavaScript

   function mandelIter(cx, cy, maxIter) {
     var x = 0.0;
     var y = 0.0;
     var xx = 0;
     var yy = 0;
     var xy = 0;
    
     var i = maxIter;
     while (i-- && xx + yy <= 4) {
       xy = x * y;
       xx = x * x;
       yy = y * y;
       x = xx - yy + cx;
       y = xy + xy + cy;
     }
     return maxIter - i;
   }
    
   function mandelbrot(canvas, xmin, xmax, ymin, ymax, iterations) {
     var width = canvas.width;
     var height = canvas.height;
    
     var ctx = canvas.getContext('2d');
     var img = ctx.getImageData(0, 0, width, height);
     var pix = img.data;
    
     for (var ix = 0; ix < width; ++ix) {
       for (var iy = 0; iy < height; ++iy) {
         var x = xmin + (xmax - xmin) * ix / (width - 1);
         var y = ymin + (ymax - ymin) * iy / (height - 1);
         var i = mandelIter(x, y, iterations);
         var ppos = 4 * (width * iy + ix);
    
         if (i > iterations) {
           pix[ppos] = 0;
           pix[ppos + 1] = 0;
           pix[ppos + 2] = 0;
         } else {
           var c = 3 * Math.log(i) / Math.log(iterations - 1.0);
    
           if (c < 1) {
             pix[ppos] = 255 * c;
             pix[ppos + 1] = 0;
             pix[ppos + 2] = 0;
           }
           else if ( c < 2 ) {
             pix[ppos] = 255;
             pix[ppos + 1] = 255 * (c - 1);
             pix[ppos + 2] = 0;
           } else {
             pix[ppos] = 255;
             pix[ppos + 1] = 255;
             pix[ppos + 2] = 255 * (c - 2);
           }
         }
         pix[ppos + 3] = 255;
       }
     }
    
     ctx.putImageData(img, 0, 0);
   }
    
   var canvas = document.createElement('canvas');
   canvas.width = 900;
   canvas.height = 600;
    
   document.body.insertBefore(canvas, document.body.childNodes[0]);
    
   mandelbrot(canvas, -2, 1, -1, 1, 1000);


Julia
=====

.. code:: Julia

   function mandelbrot(a)
       z = 0
       for i=1:50
           z = z^2 + a
       end
       return z
   end
    
   for y=1.0:-0.05:-1.0
       for x=-2.0:0.0315:0.5
           abs(mandelbrot(complex(x, y))) < 2 ? print("*") : print(" ")
       end
       println()
   end


Matlab
======

.. code:: Matlab

   function [theSet,realAxis,imaginaryAxis] = mandelbrotSet(start,gridSpacing,last,maxIteration)
    
       %Define the escape time algorithm
       function escapeTime = escapeTimeAlgorithm(z0)
    
           escapeTime = 0;
           z = 0;
    
           while( (abs(z)<=2) && (escapeTime < maxIteration) )
               z = (z + z0)^2;            
               escapeTime = escapeTime + 1;
           end
    
       end
    
       %Define the imaginary axis
       imaginaryAxis = (imag(start):imag(gridSpacing):imag(last));
    
       %Define the real axis
       realAxis = (real(start):real(gridSpacing):real(last));
    
       %Construct the complex plane from the real and imaginary axes
       complexPlane = meshgrid(realAxis,imaginaryAxis) + meshgrid(imaginaryAxis(end:-1:1),realAxis)'.*i;
    
       %Apply the escape time algorithm to each point in the complex plane 
       theSet = arrayfun(@escapeTimeAlgorithm, complexPlane);
    
    
       %Draw the set
       pcolor(realAxis,imaginaryAxis,theSet);
       shading flat;
    
   end


MySQL
=====

.. code:: MySQL

   -- Table to contain all the data points
   CREATE TABLE points (
     c_re DOUBLE,
     c_im DOUBLE,
     z_re DOUBLE DEFAULT 0,
     z_im DOUBLE DEFAULT 0,
     znew_re DOUBLE DEFAULT 0,
     znew_im DOUBLE DEFAULT 0,
     steps INT DEFAULT 0,
     active CHAR DEFAULT 1
   );
    
   DELIMITER |
    
   -- Iterate over all the points in the table 'points'
   CREATE PROCEDURE itrt (IN n INT)
   BEGIN
     label: LOOP
       UPDATE points
         SET
           znew_re=POWER(z_re,2)-POWER(z_im,2)+c_re,
           znew_im=2*z_re*z_im+c_im,
           steps=steps+1
         WHERE active=1;
       UPDATE points SET
           z_re=znew_re,
           z_im=znew_im,
           active=IF(POWER(z_re,2)+POWER(z_im,2)>4,0,1)
         WHERE active=1;
       SET n = n - 1;
       IF n > 0 THEN
         ITERATE label;
       END IF;
       LEAVE label;
     END LOOP label;
   END|
    
   -- Populate the table 'points'
   CREATE PROCEDURE populate (
     r_min DOUBLE,
     r_max DOUBLE,
     r_step DOUBLE,
     i_min DOUBLE,
     i_max DOUBLE,
     i_step DOUBLE)
   BEGIN
     DELETE FROM points;
     SET @rl = r_min;
     SET @a = 0;
     rloop: LOOP
       SET @im = i_min;
       SET @b = 0;
       iloop: LOOP
         INSERT INTO points (c_re, c_im)
           VALUES (@rl, @im);
         SET @b=@b+1;
         SET @im=i_min + @b * i_step;
         IF @im < i_max THEN
           ITERATE iloop;
         END IF;
         LEAVE iloop;
       END LOOP iloop;
         SET @a=@a+1;
       SET @rl=r_min + @a * r_step;
       IF @rl < r_max THEN
         ITERATE rloop;
       END IF;
       LEAVE rloop;
     END LOOP rloop;
   END|
    
   DELIMITER ;
    
   -- Choose size and resolution of graph
   --             R_min, R_max, R_step, I_min, I_max, I_step
   CALL populate( -2.5,  1.5,   0.005,  -2,    2,     0.005 );
    
   -- Calculate 50 iterations
   CALL itrt( 50 );
    
   -- Create the image (/tmp/image.ppm)
   -- Note, MySQL will not over-write an existing file and you may need
   -- administrator access to delete or move it
   SELECT @xmax:=COUNT(c_re) INTO @xmax FROM points GROUP BY c_im LIMIT 1;
   SELECT @ymax:=COUNT(c_im) INTO @ymax FROM points GROUP BY c_re LIMIT 1;
   SET group_concat_max_len=11*@xmax*@ymax;
   SELECT
     'P3', @xmax, @ymax, 200,
     GROUP_CONCAT(
       CONCAT(
         IF( active=1, 0, 55+MOD(steps, 200) ), ' ',
         IF( active=1, 0, 55+MOD(POWER(steps,3), 200) ), ' ',
         IF( active=1, 0, 55+MOD(POWER(steps,2), 200) ) )
       ORDER BY c_im ASC, c_re ASC SEPARATOR ' ' )
       INTO OUTFILE '/tmp/image.ppm'
     FROM points;
    


Octave
======

.. code:: Octave

   #! /usr/bin/octave -qf
   global width = 200;
   global height = 200;
   maxiter = 100;
    
   z0 = 0;
   global cmax = 1 + i;
   global cmin = -2 - i;
    
   function cs = pscale(c)
     global cmax;
     global cmin;
     global width;
     global height;
     persistent px = (real(cmax-cmin))/width;
     persistent py = (imag(cmax-cmin))/height;
     cs = real(cmin) + px*real(c) + i*(imag(cmin) + py*imag(c));
   endfunction
    
   ms = zeros(width, height);
   for x = 0:width-1
     for y = 0:height-1
       z0 = 0;
       c = pscale(x+y*i);
       for ic = 1:maxiter
         z1 = z0^2 + c;
         if ( abs(z1) > 2 ) break; endif
         z0 = z1;
       endfor
       ms(x+1, y+1) = ic/maxiter;
     endfor
   endfor
    
   saveimage("mandel.ppm", round(ms .* 255).', "ppm");


Perl
====

.. code:: Perl

   use Math::Complex;
    
   sub mandelbrot {
       my ($z, $c) = @_[0,0];
       for (1 .. 20) {
           $z = $z * $z + $c;
           return $_ if abs $z > 2;
       }
   }
    
   for (my $y = 1; $y >= -1; $y -= 0.05) {
       for (my $x = -2; $x <= 0.5; $x += 0.0315)
           {print mandelbrot($x + $y * i) ? ' ' : '#'}
       print "\n"
   }


PostScript
==========

.. code:: PostScript

   %!PS-Adobe-2.0
   %%BoundingBox: 0 0 300 200
   %%EndComments
   /origstate save def
   /ld {load def} bind def
   /m /moveto ld /g /setgray ld
   /dot { currentpoint 1 0 360 arc fill } bind def
   %%EndProlog
   % param
   /maxiter 200 def
   % complex manipulation
   /complex { 2 array astore } def
   /real { 0 get } def
   /imag { 1 get } def
   /cmul { /a exch def /b exch def
       a real b real mul
       a imag b imag mul sub
       a real b imag mul
       a imag b real mul add
       2 array astore
   } def
   /cadd { aload pop 3 -1 roll aload pop
       3 -1 roll add
       3 1 roll add exch 2 array astore
   } def
   /cconj { aload pop neg 2 array astore } def
   /cabs2 { dup cconj cmul 0 get} def
   % mandel
   200 100 translate
   -200 1 100 { /x exch def
     -100 1 100 { /y exch def
       /z0 0.0 0.0 complex def
       0 1 maxiter { /iter exch def
   	x 100 div y 100 div complex
   	z0 z0 cmul
   	cadd dup /z0 exch def
   	cabs2 4 gt {exit} if
       } for
       iter maxiter div g
       x y m dot
     } for
   } for
   %
   showpage
   origstate restore
   %%EOF


PowerShell
==========

.. code:: PowerShell

   $x = $y = $i = $j = $r = -16
   $colors = [Enum]::GetValues([System.ConsoleColor])
    
   while(($y++) -lt 15)
   {
       for($x=0; ($x++) -lt 84; Write-Host " " -BackgroundColor ($colors[$k -band 15]) -NoNewline)
       {
           $i = $k = $r = 0
    
           do
           {
               $j = $r * $r - $i * $i -2 + $x / 25
               $i = 2 * $r * $i + $y / 10
               $r = $j
           }
           while (($j * $j + $i * $i) -lt 11 -band ($k++) -lt 111)
       }
    
       Write-Host
   }


Prolog
======

.. code:: Prolog

   :- use_module(library(pce)).
    
   mandelbrot :-
       new(D, window('Mandelbrot Set')),
       send(D, size, size(700, 650)),
       new(Img, image(@nil, width := 700, height := 650, kind := pixmap)),
    
       forall(between(0,699, I),
              (   forall(between(0,649, J),
                 (   get_RGB(I, J, R, G, B),
                     R1 is (R * 256) mod 65536,
                     G1 is (G * 256) mod 65536,
                     B1 is (B * 256) mod 65536,
                     send(Img, pixel(I, J, colour(@default, R1, G1, B1))))))),
       new(Bmp, bitmap(Img)),
       send(D, display, Bmp, point(0,0)),
       send(D, open).
    
   get_RGB(X, Y, R, G, B) :-
       CX is (X - 350) / 150,
       CY is (Y - 325) / 150,
       Iter = 570,
       compute_RGB(CX, CY, 0, 0, Iter, It),
       IterF is It \/ It << 15,
       R is IterF >> 16,
       Iter1 is IterF - R << 16,
       G is Iter1 >> 8,
       B  is Iter1 - G << 8.
    
   compute_RGB(CX, CY, ZX, ZY, Iter, IterF) :-
       ZX * ZX + ZY * ZY < 4,
       Iter > 0,
       !,
       Tmp is  ZX * ZX - ZY * ZY + CX,
       ZY1 is 2 * ZX * ZY + CY,
       Iter1 is Iter - 1,
       compute_RGB(CX, CY, Tmp, ZY1, Iter1, IterF).
    
   compute_RGB(_CX, _CY, _ZX, _ZY, Iter, Iter).


Python
======

.. code:: Python

   import math
    
   def mandelbrot(z , c , n=40):
       if abs(z) > 1000:
           return float("nan")
       elif n > 0:
           return mandelbrot(z ** 2 + c, c, n - 1) 
       else:
           return z ** 2 + c
    
   print("\n".join(["".join(["#" if not math.isnan(mandelbrot(0, x + 1j * y).real) else " "
                    for x in [a * 0.02 for a in range(-80, 30)]]) 
                    for y in [a * 0.05 for a in range(-20, 20)]])
        )
    

R
=

.. code:: R

   iterate.until.escape <- function(z, c, trans, cond, max=50, response=dwell) {
     #we iterate all active points in the same array operation,
     #and keeping track of which points are still iterating.
     active <- seq_along(z)
     dwell <- z
     dwell[] <- 0
     for (i in 1:max) {
       z[active] <- trans(z[active], c[active]);
       survived <- cond(z[active])
       dwell[active[!survived]] <- i
       active <- active[survived]
       if (length(active) == 0) break
     }
     eval(substitute(response))
   }
    
   re = seq(-2, 1, len=500)
   im = seq(-1.5, 1.5, len=500)
   c <- outer(re, im, function(x,y) complex(real=x, imaginary=y))
   x <- iterate.until.escape(array(0, dim(c)), c,
                             function(z,c)z^2+c, function(z)abs(z) <= 2,
                             max=100)
   image(x)


Ruby
====

.. code:: Ruby

   require 'complex'
    
   def mandelbrot(a)
     Array.new(50).inject(0) { |z,c| z*z + a }
   end
    
   (1.0).step(-1,-0.05) do |y|
     (-2.0).step(0.5,0.0315) do |x|
       print mandelbrot(Complex(x,y)).abs < 2 ? '*' : ' '
     end
     puts
   end


Scala
=====

.. code:: Scala

   import org.rosettacode.ArithmeticComplex._
   import java.awt.Color
    
   object Mandelbrot
   {
      def generate(width:Int =600, height:Int =400)={
         val bm=new RgbBitmap(width, height)
    
         val maxIter=1000
         val xMin = -2.0
         val xMax =  1.0
         val yMin = -1.0
         val yMax =  1.0
    
         val cx=(xMax-xMin)/width
         val cy=(yMax-yMin)/height
    
         for(y <- 0 until bm.height; x <- 0 until bm.width){
            val c=Complex(xMin+x*cx, yMin+y*cy)
            val iter=itMandel(c, maxIter, 4)
            bm.setPixel(x, y, getColor(iter, maxIter))
         }
         bm
      }
    
      def itMandel(c:Complex, imax:Int, bailout:Int):Int={
         var z=Complex()
         for(i <- 0 until imax){
            z=z*z+c;
            if(z.abs > bailout) return i
         }
         imax;
      }
    
      def getColor(iter:Int, max:Int):Color={
         if (iter==max) return Color.BLACK
    
         var c=3*math.log(iter)/math.log(max-1.0)
         if(c<1) new Color((255*c).toInt, 0, 0)
         else if(c<2) new Color(255, (255*(c-1)).toInt, 0)
         else new Color(255, 255, (255*(c-2)).toInt)
      }
   }


Scheme
======

.. code:: Scheme

   (define x-centre -0.5)
   (define y-centre 0.0)
   (define width 4.0)
   (define i-max 800)
   (define j-max 600)
   (define n 100)
   (define r-max 2.0)
   (define file "out.pgm")
   (define colour-max 255)
   (define pixel-size (/ width i-max))
   (define x-offset (- x-centre (* 0.5 pixel-size (+ i-max 1))))
   (define y-offset (+ y-centre (* 0.5 pixel-size (+ j-max 1))))
    
   (define (inside? z)
     (define (*inside? z-0 z n)
       (and (< (magnitude z) r-max)
            (or (= n 0)
                (*inside? z-0 (+ (* z z) z-0) (- n 1)))))
     (*inside? z 0 n))
    
   (define (boolean->integer b)
     (if b colour-max 0))
    
   (define (pixel i j)
     (boolean->integer
       (inside?
         (make-rectangular (+ x-offset (* pixel-size i))
                           (- y-offset (* pixel-size j))))))
    
   (define (plot)
     (with-output-to-file file
       (lambda ()
         (begin (display "P2") (newline)
                (display i-max) (newline)
                (display j-max) (newline)
                (display colour-max) (newline)
                (do ((j 1 (+ j 1))) ((> j j-max))
                    (do ((i 1 (+ i 1))) ((> i i-max))
                        (begin (display (pixel i j)) (newline))))))))
    
   (plot)


Tcl
===

.. code:: Tcl

   package require Tk
    
   proc mandelIters {cx cy} {
       set x [set y 0.0]
       for {set count 0} {hypot($x,$y) < 2 && $count < 255} {incr count} {
           set x1 [expr {$x*$x - $y*$y + $cx}]
           set y1 [expr {2*$x*$y + $cy}]
           set x $x1; set y $y1
       }
       return $count
   }
   proc mandelColor {iter} {
       set r [expr {16*($iter % 15)}]
       set g [expr {32*($iter % 7)}]
       set b [expr {8*($iter % 31)}]
       format "#%02x%02x%02x" $r $g $b
   }
   image create photo mandel -width 300 -height 300
   # Build picture in strips, updating as we go so we have "progress" monitoring
   # Also set the cursor to tell the user to wait while we work.
   pack [label .mandel -image mandel -cursor watch]
   update
   for {set x 0} {$x < 300} {incr x} {
       for {set y 0} {$y < 300} {incr y} {
           set i [mandelIters [expr {($x-220)/100.}] [expr {($y-150)/90.}]]
           mandel put [mandelColor $i] -to $x $y
       }
       update
   }
   .mandel configure -cursor {}


XSLT
====

.. code:: XSLT

   <?xml version="1.0" encoding="UTF-8"?>
   <xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">
    
   <!-- XSLT Mandelbrot - written by Joel Yliluoma 2007, http://iki.fi/bisqwit/ -->
    
   <xsl:output method="html" indent="no"
     doctype-public="-//W3C//DTD HTML 4.01//EN"
     doctype-system="http://www.w3.org/TR/REC-html40/strict.dtd"
    />
    
   <xsl:template match="/fractal">
    <html>
     <head>
      <title>XSLT fractal</title>
      <style type="text/css">
   body { color:#55F; background:#000 }
   pre { font-family:monospace; font-size:7px }
   pre span { background:<xsl:value-of select="background" /> }
      </style>
     </head>
     <body>
      <div style="position:absolute;top:20px;left:20em">
       Copyright © 1992,2007 Joel Yliluoma
       (<a href="http://iki.fi/bisqwit/">http://iki.fi/bisqwit/</a>)
      </div>
      <h1 style="margin:0px">XSLT fractal</h1>
      <pre><xsl:call-template name="bisqwit-mandelbrot" /></pre>
     </body>
    </html>
   </xsl:template>
    
   <xsl:template name="bisqwit-mandelbrot"
     ><xsl:call-template name="bisqwit-mandelbrot-line">
      <xsl:with-param name="y" select="y/min"/>
     </xsl:call-template
   ></xsl:template>
    
   <xsl:template name="bisqwit-mandelbrot-line"
    ><xsl:param name="y"
    /><xsl:call-template name="bisqwit-mandelbrot-column">
     <xsl:with-param name="x" select="x/min"/>
     <xsl:with-param name="y" select="$y"/>
    </xsl:call-template
    ><xsl:if test="$y < y/max"
     ><br
     /><xsl:call-template name="bisqwit-mandelbrot-line">
      <xsl:with-param name="y" select="$y + y/step"/>
     </xsl:call-template
    ></xsl:if
   ></xsl:template>
    
   <xsl:template name="bisqwit-mandelbrot-column"
    ><xsl:param name="x"
    /><xsl:param name="y"
    /><xsl:call-template name="bisqwit-mandelbrot-slot">
     <xsl:with-param name="x" select="$x" />
     <xsl:with-param name="y" select="$y" />
     <xsl:with-param name="zr" select="$x" />
     <xsl:with-param name="zi" select="$y" />
    </xsl:call-template
    ><xsl:if test="$x < x/max"
     ><xsl:call-template name="bisqwit-mandelbrot-column">
      <xsl:with-param name="x" select="$x + x/step"/>
      <xsl:with-param name="y" select="$y" />
     </xsl:call-template
    ></xsl:if
   ></xsl:template>
    
   <xsl:template name="bisqwit-mandelbrot-slot"
   ><xsl:param name="x"
    /><xsl:param name="y"
    /><xsl:param name="zr"
    /><xsl:param name="zi"
    /><xsl:param name="iter" select="0"
    /><xsl:variable name="zrsqr" select="($zr * $zr)"
    /><xsl:variable name="zisqr" select="($zi * $zi)"
    /><xsl:choose>
     <xsl:when test="(4*scale*scale >= $zrsqr + $zisqr) and (maxiter > $iter+1)"
      ><xsl:call-template name="bisqwit-mandelbrot-slot">
       <xsl:with-param name="x" select="$x" />
       <xsl:with-param name="y" select="$y" />
       <xsl:with-param name="zi" select="(2 * $zr * $zi) div scale + $y" />
       <xsl:with-param name="zr" select="($zrsqr - $zisqr) div scale + $x" />
       <xsl:with-param name="iter" select="$iter + 1" />
      </xsl:call-template
     ></xsl:when>
     <xsl:otherwise
      ><xsl:variable name="magnitude" select="magnitude[@value=$iter]"
       /><span style="color:{$magnitude/color}"
      ><xsl:value-of select="$magnitude/symbol"
     /></span></xsl:otherwise>
    </xsl:choose
   ></xsl:template>
    
   </xsl:stylesheet>
    

End Thought
===========

There is more than one way to accomplish almost anything.

