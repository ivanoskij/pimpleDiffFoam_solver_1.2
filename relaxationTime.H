//////////////////////////////////////////////////////////////////////////////////////

Info<< "Reading Schmidt number\n" << endl;
dimensionedScalar Sc
(
    "Sc",
    transportProperties.lookup("Sc")
);
//////////////////////////////////////////////////////////////////////////////////////

Info<< "Reading droplet diameter\n" << endl;
dimensionedScalar dW
(
    "dW",
    transportProperties.lookup("dW")
);
Info<< "Reading water density\n" << endl;
dimensionedScalar rhoW
(
    "rhoW",
    transportProperties.lookup("rhoW")
);
Info<< "Reading air viscosity\n" << endl;
dimensionedScalar mu
(
    "mu",
    transportProperties.lookup("mu")
);
Info<< "Reading air kinematic viscosity\n" << endl;
dimensionedScalar nu
(
    "nu",
    transportProperties.lookup("nu")
);

//////////////////////////////////////////////////////////////////////////////////////

const dimensionedScalar dummyV("dummyV",dimensionSet(0, 1, -1, 0, 0),1e-10);

Info<< "Computing Reloc\n" << endl;
volScalarField Reloc ("Reloc", ((max(mag(dummyV),mag(U)))*dW*rhoW/mu));		//local reynolds number//[adim]

Info<< "Computing Cw\n" << endl;
volScalarField Cw ("Cw", (24/Reloc+4*pow(Reloc,-0.333)));			//drag coefficient//[adim]

Info<< "Comuting relaxation time\n" << endl;
volScalarField relt ("relt",(Cw*rhoW*pow(dW,2)/(18*mu)));			//relaxation time//[s]

//////////////////////////////////////////////////////////////////////////////////////

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

const dimensionedVector g(gravitationalProperties.lookup("g"));
const dimensionedScalar magg = mag(g);
const dimensionedVector gHat = g/magg;

Info<< "Reading field grav\n" << endl;
volVectorField grav
(
    IOobject
    (
        "grav",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    g
);

//////////////////////////////////////////////////////////////////////////////////////

Info<< "Computing settling velocity\n" << endl;
volVectorField Uw
(
    IOobject
    (
        "Uw",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    relt*grav
);

//////////////////////////////////////////////////////////////////////////////////////

Info<< "computing face flux field phiW\n" << endl;
surfaceScalarField phiW
(
    IOobject
    (
        "phiW",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    fvc::flux(Uw)
);
//////////////////////////////////////////////////////////////////////////////////////
