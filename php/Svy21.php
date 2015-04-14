<?php

class Svy21{
    var $a = 6378137;
    var $f;

    // SVY21 Projection
    // Fundamental point: Base 7 at Pierce Resevoir.
    // Latitude: 1 22 02.9154 N, longitude: 103 49 31.9752 E (of Greenwich).

    // Known Issue: Setting (oLat, oLon) to the exact coordinates specified above
		// results in computation being slightly off. The values below give the most 
    // accurate represenation of test data.
    var $oLat = 1.366666;     // origin's lat in degrees
    var $oLon = 103.833333;   // origin's lon in degrees
    var $oN = 38744.572;      // false Northing
	var $oE = 28001.642;      // false Easting
    var $k = 1;               // scale factor


    public function computeSVY21($lat, $lon){
        //Returns a pair (N, E) representing Northings and Eastings in SVY21.

         $latR = $lat * pi() / 180;
         $sinLat = sin($latR);
         $sin2Lat = $sinLat * $sinLat;
         $cosLat = cos($latR);
         $cos2Lat = $cosLat * $cosLat;
         $cos3Lat = $cos2Lat * $cosLat;
         $cos4Lat = $cos3Lat * $cosLat;
         $cos5Lat = $cos4Lat * $cosLat;
         $cos6Lat = $cos5Lat * $cosLat;
         $cos7Lat = $cos6Lat * $cosLat;

         $rho = $this->calcRho($sin2Lat);
         $v = $this->calcV($sin2Lat);
         $psi = $v / $rho;
         $t = tan($latR);
         $w = ($lon - $this->oLon) * pi() / 180;

         $M = $this->calcM($lat);
         $Mo = $this->calcM($this->oLat);

         $w2 = $w * $w;
         $w4 = $w2 * $w2;
         $w6 = $w4 * $w2;
         $w8 = $w6 * $w2;

         $psi2 = $psi * $psi;
         $psi3 = $psi2 * $psi;
         $psi4 = $psi3 * $psi;

         $t2 = $t * $t;
         $t4 = $t2 * $t2;
         $t6 = $t4 * $t2;

        //	Compute Northing
         $nTerm1 = $w2 / 2 * $v * $sinLat * $cosLat;
         $nTerm2 = $w4 / 24 * $v * $sinLat * $cos3Lat * (4 * $psi2 + $psi - $t2);
         $nTerm3 = $w6 / 720 * $v * $sinLat * $cos5Lat * ((8 * $psi4) * (11 - 24 * $t2) - (28 * $psi3) * (1 - 6 * $t2) + $psi2 * (1 - 32 * $t2) - $psi * 2 * $t2 + $t4);
         $nTerm4 = $w8 / 40320 * $v * $sinLat * $cos7Lat * (1385 - 3111 * $t2 + 543 * $t4 - $t6);
         $N = $this->oN + $this->k * ($M - $Mo + $nTerm1 + $nTerm2 + $nTerm3 + $nTerm4);

        //	Compute Easting
         $eTerm1 = $w2 / 6 * $cos2Lat * ($psi - $t2);
         $eTerm2 = $w4 / 120 * $cos4Lat * ((4 * $psi3) * (1 - 6 * $t2) + $psi2 * (1 + 8 * $t2) - $psi * 2 * $t2 + $t4);
         $eTerm3 = $w6 / 5040 * $cos6Lat * (61 - 479 * $t2 + 179 * $t4 - $t6);
         $E = $this->oE + $this->k * $v * $w * $cosLat * (1 + $eTerm1 + $eTerm2 + $eTerm3);

        return [$E,$N];
	}

		
		
	private function calcM($lat){
         $latR = $lat * pi() / 180;
        return $this->a * (($this->A0 * $latR) - ($this->A2 * sin(2 * $latR)) + ($this->A4 * sin(4 * $latR)) - ($this->A6 * sin(6 * $latR)));
		}
				
    private function calcRho($sin2Lat){
         $num = $this->a * (1 - $this->e2);
         $denom = pow(1 - $this->e2 * $sin2Lat, 3. / 2.);
        return $num / $denom;
		}

    private function calcV($sin2Lat){
        $poly = 1 - $this->e2 * $sin2Lat;
        return $this->a / sqrt($poly);
		}

	public function __construct()
	{
		$this->f = 1 / 298.257223563;
		$this->b = $this->a * (1 - $this->f);
        $this->e2 = (2 * $this->f) - ($this->f * $this->f);
        $this->e4 = $this->e2 * $this->e2;
        $this->e6 = $this->e4 * $this->e2;
        $this->A0 = 1 - ($this->e2 / 4) - (3 * $this->e4 / 64) - (5 * $this->e6 / 256);
        $this->A2 = (3. / 8.) * ($this->e2 + ($this->e4 / 4) + (15 * $this->e6 / 128));
        $this->A4 = (15. / 256.) * ($this->e4 + (3 * $this->e6 / 4));
        $this->A6 = 35 * $this->e6 / 3072;
		
	}

}
