stk.v.12.0
WrittenBy    STK_v12.2.0

BEGIN Planet

    Name		 Jupiter

    BEGIN PathDescription

        CentralBody		 Jupiter
        UseCbEphemeris		 Yes

        BEGIN EphemerisData

            EphemerisSource		 None

            JplIndex		 4

            JplSpiceId		 599

            ApplyTDTtoTDBCorrectionForDE		 Yes

            OrbitEpoch		  2.4515450000000000e+06
            OrbitCrdSys		 EclipticJ2000ICRF
            OrbitCrdSysEpoch		  2.4515450000000000e+06
            OrbitParentGM		  1.3271244004193939e+20
            OrbitMeanDist		  7.7834081669271082e+11
            OrbitEcc		  4.8386239999999997e-02
            OrbitInc		  1.3043969500000001e+00
            OrbitRAAN		  1.0047390909000001e+02
            OrbitPerLong		  1.4728479829999999e+01
            OrbitMeanLong		  3.4396440510000005e+01
            OrbitMeanDistDot		 -4.7539561539080086e+02
            OrbitEccDot		 -3.6284736481861740e-09
            OrbitIncDot		 -5.0298151950718684e-08
            OrbitRAANDot		  5.6041357973990419e-06
            OrbitPerLongDot		  5.8186633812457225e-06
            OrbitMeanLongDot		  8.3086820746064355e-02

        END EphemerisData

    END PathDescription

    BEGIN PhysicalData

        GM		  1.2668653500000000e+17
        Radius		  7.1492000000000000e+07
        Magnitude		  0.0000000000000000e+00
        ReferenceDistance		  0.0000000000000000e+00

    END PhysicalData

    BEGIN AutoRename

        AutoRename		 Yes

    END AutoRename

    BEGIN Extensions

        BEGIN ExternData
        END ExternData

        BEGIN ADFFileData
        END ADFFileData

        BEGIN AccessConstraints
            LineOfSight IncludeIntervals

            UsePreferredMaxStep No
            PreferredMaxStep 360
        END AccessConstraints

        BEGIN Desc
            BEGIN ShortText

            END ShortText
            BEGIN LongText

            END LongText
        END Desc

        BEGIN Crdn
        END Crdn

        BEGIN Graphics

            BEGIN Attributes

                MarkerColor		 #ffff00
                LabelColor		 #ffff00
                LineColor		 #ffff00
                LineStyle		 0
                LineWidth		 1
                MarkerStyle		 2
                FontStyle		 0

            END Attributes

            BEGIN Graphics

                Show		 On
                Inherit		 On
                ShowLabel		 On
                ShowPlanetPoint		 On
                ShowSubPlanetPoint		 Off
                ShowSubPlanetLabel		 Off
                ShowOrbit		 On
                NumOrbitPoints		 360
                OrbitTime		  0.0000000000000000e+00
                OrbitDisplay		                OneOrbit		
                TransformTrajectory		 On

            END Graphics
        END Graphics

        BEGIN VO
        END VO

    END Extensions

END Planet

