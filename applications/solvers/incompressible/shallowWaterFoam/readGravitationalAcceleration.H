Info<< "\nReading gravitationalProperties" << endl;

IOdictionary gravitationalProperties
(
    IOobject
    (
        "gravitationalProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

const dimensionedVector g("g", dimAcceleration, gravitationalProperties);
const bool rotating(gravitationalProperties.get<bool>("rotating"));
const dimensionedVector Omega =
(
    rotating
  ? dimensionedVector("Omega", dimless/dimTime, gravitationalProperties)
  : dimensionedVector("Omega", dimless/dimTime, Zero)
);
const dimensionedScalar magg = mag(g);
const dimensionedVector gHat = g/magg;
