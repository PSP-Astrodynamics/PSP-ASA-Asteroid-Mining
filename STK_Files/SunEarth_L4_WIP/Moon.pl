stk.v.12.0
WrittenBy    STK_v12.2.0

BEGIN Planet

    Name		 Moon

    BEGIN PathDescription

        CentralBody		 Moon
        UseCbEphemeris		 Yes

        BEGIN EphemerisData

            EphemerisSource		 JplDE

            JplIndex		 9

            JplSpiceId		 301

            ApplyTDTtoTDBCorrectionForDE		 Yes

            OrbitEpoch		  2.4515450000000000e+06
            OrbitCrdSys		 ICRF
            OrbitCrdSysEpoch		  2.4515450000000000e+06
            OrbitParentGM		  1.3271244004193939e+20
            OrbitMeanDist		  1.4959826115044250e+11
            OrbitEcc		  1.6711230000000001e-02
            OrbitInc		 -1.5310000000000001e-05
            OrbitRAAN		  0.0000000000000000e+00
            OrbitPerLong		  1.0293768193000000e+02
            OrbitMeanLong		  1.0046457166000002e+02
            OrbitMeanDistDot		  2.3018207620369612e+01
            OrbitEccDot		 -1.2024640657084188e-09
            OrbitIncDot		 -3.5446078028747444e-07
            OrbitRAANDot		  0.0000000000000000e+00
            OrbitPerLongDot		  8.8507498973305977e-06
            OrbitMeanLongDot		  9.8560910197973994e-01

        END EphemerisData

    END PathDescription

    BEGIN PhysicalData

        GM		  4.9028003055554004e+12
        Radius		  1.7374000000000000e+06
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
        END Desc

        BEGIN Crdn
        END Crdn

        BEGIN Graphics

            BEGIN Attributes

                MarkerColor		 #0eff82
                LabelColor		 #0eff82
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
                OrbitTime		  2.3806552530814158e+06
                OrbitDisplay		                OneOrbit		
                TransformTrajectory		 On

            END Graphics
        END Graphics

        BEGIN VO
        END VO

    END Extensions

END Planet

