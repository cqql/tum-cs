within ;
package Hydro "Hydraulic Components"

  import SI = Modelica.SIunits;

  package Examples

    model SourceExample

      Sources.ConstantPressure pressureSource(p=1000000)
        annotation (Placement(transformation(extent={{-52,6},{-32,26}})));
      Basic.Pipe pipe(k=1)
        annotation (Placement(transformation(extent={{-2,6},{18,26}})));
      Basic.CompressibleVolume compressibleVolume(k=1, V=2)
        annotation (Placement(transformation(extent={{46,6},{66,26}})));
    equation
      connect(pressureSource.out, pipe.inFlow) annotation (Line(
          points={{-35.4,16},{0,16},{0,16}},
          color={0,127,0},
          thickness=1));
      connect(pipe.outFlow, compressibleVolume.inout) annotation (Line(
          points={{16,16},{34,16},{50,16}},
          color={0,127,0},
          thickness=1));
      annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
                -100},{100,100}})));
    end SourceExample;

    model SpeedOfSound

      Sources.SignalPressure signalPressure
        annotation (Placement(transformation(extent={{-74,-10},{-54,10}})));
      AirTube airTube
        annotation (Placement(transformation(extent={{-44,-10},{-24,10}})));
      AirTube airTube1
        annotation (Placement(transformation(extent={{-14,-10},{6,10}})));
      AirTube airTube2
        annotation (Placement(transformation(extent={{8,-10},{28,10}})));
      AirTube airTube3
        annotation (Placement(transformation(extent={{36,-10},{56,10}})));
      AirTube airTube4
        annotation (Placement(transformation(extent={{68,-10},{88,10}})));
      Modelica.Blocks.Sources.Sine sine(amplitude=100, freqHz=60)
        annotation (Placement(transformation(extent={{-110,-10},{-90,10}})));
    equation
      connect(airTube.positivePin, airTube1.negativePin) annotation (Line(
          points={{-26,0},{-26,0},{-12,0}},
          color={0,127,0},
          thickness=1));
      connect(airTube1.positivePin, airTube2.negativePin) annotation (Line(
          points={{4,0},{7,0},{10,0}},
          color={0,127,0},
          thickness=1));
      connect(airTube2.positivePin, airTube3.negativePin) annotation (Line(
          points={{26,0},{32,0},{38,0}},
          color={0,127,0},
          thickness=1));
      connect(airTube3.positivePin, airTube4.negativePin) annotation (Line(
          points={{54,0},{62,0},{70,0}},
          color={0,127,0},
          thickness=1));
      connect(sine.y, signalPressure.signal) annotation (Line(points={{-89,0},{
              -78,0},{-78,0},{-70,0}}, color={0,0,127}));
      connect(signalPressure.outFlow, airTube.negativePin) annotation (Line(
          points={{-56,0},{-42,0},{-42,0}},
          color={0,127,0},
          thickness=1));
      annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{
                -100,-100},{100,100}})));
    end SpeedOfSound;

    model AirTube

      Interfaces.NegativePin negativePin
        annotation (Placement(transformation(extent={{-90,-10},{-70,10}})));
      Interfaces.PositivePin positivePin
        annotation (Placement(transformation(extent={{70,-10},{90,10}})));
      Basic.Inductor inductor(
        l=0.2,
        rho=1.29,
        A=0.01)
        annotation (Placement(transformation(extent={{-48,-10},{-28,10}})));
      Basic.Pipe pipe(k=5000)
        annotation (Placement(transformation(extent={{10,-10},{30,10}})));
      Basic.CompressibleVolume compressibleVolume(k=0.00001, V=0.002)
        annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=90,
            origin={-10,20})));
    equation
      connect(inductor.outFlow, pipe.inFlow) annotation (Line(
          points={{-30,0},{-30,0},{12,0}},
          color={0,127,0},
          thickness=1));
      connect(compressibleVolume.inout, pipe.inFlow) annotation (Line(
          points={{-10,14},{-10,0},{12,0}},
          color={170,255,85},
          thickness=1));
      connect(pipe.outFlow, positivePin) annotation (Line(
          points={{28,0},{80,0}},
          color={0,127,0},
          thickness=1));
      connect(negativePin, inductor.inFlow) annotation (Line(
          points={{-80,0},{-46,0}},
          color={170,255,85},
          thickness=1));
      annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{
                -100,-100},{100,100}})));
    end AirTube;

    model SignalPressureSineInteraction

      Sources.SignalPressure signalPressure
        annotation (Placement(transformation(extent={{8,-10},{28,10}})));
      Modelica.Blocks.Sources.Sine sine(amplitude=3, freqHz=2)
        annotation (Placement(transformation(extent={{-52,-10},{-32,10}})));
    equation
      connect(sine.y, signalPressure.signal)
        annotation (Line(points={{-31,0},{-10,0},{12,0}}, color={0,0,127}));
      annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{
                -100,-100},{100,100}})));
    end SignalPressureSineInteraction;
  end Examples;

  package Basic

    model Pipe

      parameter Real k "Dissipation rate";

      Interfaces.NegativePin inFlow
        annotation (Placement(transformation(extent={{-90,-10},{-70,10}})));
      Interfaces.PositivePin outFlow
        annotation (Placement(transformation(extent={{70,-10},{90,10}})));

    equation
      inFlow.p - outFlow.p = k * inFlow.v;
      outFlow.v + inFlow.v = 0;

      annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
                -100},{100,100}})), Icon(graphics={
                                             Rectangle(
              extent={{-66,26},{66,-26}},
              lineThickness=1,
              fillPattern=FillPattern.Solid,
              fillColor={170,213,255},
              pattern=LinePattern.None),
            Line(
              points={{-66,40},{66,40}},
              color={95,95,95},
              thickness=0.5),
            Line(
              points={{-66,-40},{66,-40}},
              color={95,95,95},
              thickness=0.5)}));
    end Pipe;

    model CompressibleVolume

      Interfaces.NegativePin inout annotation (Placement(transformation(extent={{-70,
                -10},{-50,10}}), iconTransformation(extent={{-70,-10},{-50,10}})));

      parameter Real k "Compressibility";
      parameter SI.Volume V "Volume";

    equation
      der(inout.p) * k * V = inout.v;

      annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},
                {100,100}}), graphics={Line(
              points={{-60,20},{-20,20},{-20,60},{60,60},{60,-60},{-20,-60},{-20,-20},
                  {-60,-20}},
              color={255,170,85},
              thickness=1), Rectangle(
              extent={{0,40},{40,-40}},
              lineThickness=1,
              fillColor={170,213,255},
              fillPattern=FillPattern.Solid,
              pattern=LinePattern.None)}));
    end CompressibleVolume;

    model Inductor

      parameter SI.Length l "Length of the pipe";
      parameter SI.Area A "Inlet/outlet area";
      parameter Real rho "Volumetric density";

      Interfaces.NegativePin inFlow
        annotation (Placement(transformation(extent={{-90,-10},{-70,10}})));
      Interfaces.PositivePin outFlow
        annotation (Placement(transformation(extent={{70,-10},{90,10}})));
      annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
                -100},{100,100}})), Icon(coordinateSystem(preserveAspectRatio=false,
              extent={{-100,-100},{100,100}}),
                                         graphics={
            Line(
              points={{-60,60},{60,60}},
              color={95,95,95},
              thickness=1),
            Rectangle(
              extent={{-60,40},{60,-40}},
              pattern=LinePattern.None,
              lineThickness=1,
              fillColor={0,128,255},
              fillPattern=FillPattern.Solid),
            Line(
              points={{-60,-60},{60,-60}},
              color={95,95,95},
              thickness=1),
            Ellipse(
              extent={{-34,24},{36,-22}},
              pattern=LinePattern.None,
              lineThickness=0.5,
              fillColor={255,255,0},
              fillPattern=FillPattern.Solid)}));

    equation
      inFlow.p - outFlow.p = der(inFlow.v) * rho * l / A;
      inFlow.v + outFlow.v = 0;

    end Inductor;
  end Basic;

  package Interfaces

    connector PositivePin
      SI.Pressure p;
      flow SI.VolumeFlowRate v;
      annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},
                {100,100}}), graphics={Rectangle(
              extent={{-80,80},{80,-80}},
              lineColor={0,127,0},
              lineThickness=1,
              fillColor={170,255,85},
              fillPattern=FillPattern.Solid), Text(
              extent={{-80,80},{80,-80}},
              lineColor={255,255,255},
              lineThickness=1,
              fillPattern=FillPattern.Solid,
              textString="+",
              textStyle={TextStyle.Bold})}));
    end PositivePin;

    connector NegativePin
      SI.Pressure p;
      flow SI.VolumeFlowRate v;
      annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,
                -100},{100,100}}), graphics={Rectangle(
              extent={{-80,80},{80,-80}},
              lineColor={170,255,85},
              fillColor={0,127,0},
              fillPattern=FillPattern.Solid,
              lineThickness=1), Text(
              extent={{80,-80},{-80,80}},
              lineColor={255,255,255},
              lineThickness=1,
              fillPattern=FillPattern.Solid,
              textString="-",
              textStyle={TextStyle.Bold})}));
    end NegativePin;
  end Interfaces;

  package Sources

    model ConstantPressure

      parameter SI.Pressure p;

      Interfaces.PositivePin out
        annotation (Placement(transformation(extent={{56,-10},{76,10}})));

    equation
      out.p = p;

      annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
                -100},{100,100}})), Icon(graphics={
                                             Rectangle(
              extent={{-80,80},{60,-80}},
              lineColor={255,0,0},
              lineThickness=1,
              fillPattern=FillPattern.Solid,
              fillColor={255,0,0}), Text(
              extent={{-36,50},{26,-50}},
              lineColor={0,255,0},
              lineThickness=1,
              fillPattern=FillPattern.Solid,
              textString="P+")}));
    end ConstantPressure;

    model SignalPressure

      Interfaces.PositivePin outFlow annotation (Placement(transformation(extent={{70,
                -10},{90,10}}), iconTransformation(extent={{70,-10},{90,10}})));
      Modelica.Blocks.Interfaces.RealInput signal annotation (Placement(
            transformation(extent={{-80,-20},{-40,20}}), iconTransformation(extent={
                {-80,-20},{-40,20}})));
      annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},
                {100,100}}), graphics={
            Polygon(
              points={{60,0},{60,0}},
              pattern=LinePattern.None,
              lineThickness=1,
              fillColor={0,128,255},
              fillPattern=FillPattern.Solid),
            Polygon(
              points={{60,0},{32,18},{60,0}},
              pattern=LinePattern.None,
              lineThickness=1,
              fillColor={0,128,255},
              fillPattern=FillPattern.Solid),
            Polygon(
              points={{60,0},{-20,60},{-20,-60},{60,0}},
              pattern=LinePattern.None,
              lineThickness=1,
              fillColor={0,128,255},
              fillPattern=FillPattern.Solid),
            Line(
              points={{60,20},{-20,80}},
              color={95,95,95},
              thickness=1),
            Line(
              points={{60,-20},{-20,-80}},
              color={95,95,95},
              thickness=1)}));

    equation
      outFlow.p = signal;

    end SignalPressure;
  end Sources;
  annotation (uses(Modelica(version="3.2.1")));
end Hydro;
