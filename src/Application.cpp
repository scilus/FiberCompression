
#include<QElapsedTimer>
#include<QFile>

#include "Application.h"


/********************************************//**
\brief default constructor
\param parent : parent of the mainwindow
***********************************************/
Application::Application(int argc, char** argv) :
QCoreApplication(argc, argv),
m_pFibers(NULL),
m_pThread(NULL),
m_inputArg(NULL),
m_outputArg(NULL),
m_compArg(NULL),
m_decompArg(NULL),
m_errorArg(NULL),
m_ttypeArg(NULL),
m_qtypeArg(NULL),
m_qnumArg(NULL),
m_etypeArg(NULL),
m_statsArg(NULL)
{
    // Init and parse commandline options
    if(initAndParseCommandLineOptions(argc, argv))
    {
    
        // Initialize Fibers instance
        m_pFibers = new Fibers();
        
        // Create thread for processing
        m_pThread = new QThread();
        m_pFibers->moveToThread(m_pThread);
        
        // Connect signals and slots
        connectSlots();
        
        // Get commandline values and validate
        if(getAndValidateCommandLineArgs(m_param))
        {
            // Validate commandline
            startFromCommandLine(m_param);
        }
        else
        {
            // Invalid argument : close application
            close();
        }
    }
    else
    {
        close();
    }
}


/********************************************//**
\brief default destructor
***********************************************/
Application::~Application()
{
    // Free pointers
    if(m_pThread != NULL)
    {
        delete m_pThread;
        m_pThread = NULL;
    }
    if(m_inputArg != NULL)
    {
        delete m_inputArg;
        m_inputArg = NULL;
    }
    if(m_outputArg != NULL)
    {
        delete m_outputArg;
        m_outputArg = NULL;
    }
    if(m_compArg != NULL)
    {
        delete m_compArg;
        m_compArg = NULL;
    }
    if(m_decompArg != NULL)
    {
        delete m_decompArg;
        m_decompArg = NULL;
    }
    if(m_errorArg != NULL)
    {
        delete m_errorArg;
        m_errorArg = NULL;
    }
    if(m_ttypeArg != NULL)
    {
        delete m_ttypeArg;
        m_ttypeArg = NULL;
    }
    if(m_qtypeArg != NULL)
    {
        delete m_qtypeArg;
        m_qtypeArg = NULL;
    }
    if(m_qnumArg != NULL)
    {
        delete m_qnumArg;
        m_qnumArg = NULL;
    }
    if(m_etypeArg != NULL)
    {
        delete m_etypeArg;
        m_etypeArg = NULL;
    }
    if(m_statsArg != NULL)
    {
        delete m_statsArg;
        m_statsArg = NULL;
    }
}


/********************************************//**
\brief connect signals and slots
***********************************************/
void Application::connectSlots()
{
    connect(m_pFibers, SIGNAL(showMessage(const QString&)), this, SLOT(showMessage(const QString&)), Qt::DirectConnection);
    connect(m_pFibers, SIGNAL(isStopped()), this, SLOT(stopProcess()), Qt::DirectConnection);
    connect(m_pFibers, SIGNAL(finished()), this, SLOT(processFinished()), Qt::DirectConnection);
    connect(m_pFibers, SIGNAL(finished()), m_pFibers, SLOT(deleteLater()));
    connect(m_pThread, SIGNAL(started()), m_pFibers, SLOT(process()));
    connect(m_pThread, SIGNAL(finished()), m_pThread, SLOT(deleteLater()));
    connect(this, SIGNAL(aboutToQuit()), this, SLOT(close()));
}


void Application::stopProcess()
{
    close();
    
}


/********************************************//**
\brief slot triggered when the status bar message must change
***********************************************/
void Application::showMessage(const QString& message)
{
    std::cout << message.toStdString() << std::endl;
}


/********************************************//**
\brief initialize commandline options and parse
\param argc : argument count (coming from main)
\param argv : argument vector (coming from main)
***********************************************/
bool Application::initAndParseCommandLineOptions(int argc, char** argv)
{
    try
    {
        
        TCLAP::CmdLine cmd("Fibercompression : \n   Tool to compress and decompress tractography files. Supported formats are .tck, .trk and .zfib.", ' ', "1.0", true);
        
        m_inputArg = new TCLAP::UnlabeledValueArg<std::string>("input",
                                                               "Input file path (supported formats are .tck, .trk and .zfib)",
                                                               true, "", "string");
        m_outputArg = new TCLAP::UnlabeledValueArg<std::string>("output",
                                                                "Output file path (supported formats are .tck, .trk and .zfib)",
                                                                true, "", "string");
        m_compArg = new TCLAP::SwitchArg("c", "compress", "Flag to compress", false);
        m_decompArg = new TCLAP::SwitchArg("d", "decompress", "Flag to decompress", false);
        
        m_forceArg = new TCLAP::SwitchArg("f", "", "Flag to overwrite existing output file", false);
        
        m_errorArg = new TCLAP::ValueArg<float>("e", "",
                                                "Maximum error in mm (default " + valToString(DEFAULT_MAX_ERROR) + ")",
                                                false, 0.5, "float");
        // Add arguments to commandline
        cmd.add(m_inputArg);
        cmd.add(m_outputArg);
        cmd.xorAdd(m_compArg, m_decompArg);
        cmd.add(m_forceArg);
        cmd.add(m_errorArg);
        

        #if !_USE_LIGHTWEIGHT
            // Create Transformation description string
            std::vector<int> transfVal;
            std::string transfStr = "Transformation type : ";
            for(TransformationType i = DCT_TRANSFORMATION; i <= NO_TRANSFORMATION; i++)
            {
                transfStr += getTransformationStringFromEnum(i).toStdString() + " : " + valToString(i);
                
                if(i == DEFAULT_TRANSF_TYPE)
                {
                    transfStr += " (default)";
                }
                if(i < NO_TRANSFORMATION)
                {
                    transfStr += ", ";
                }
                transfVal.push_back(i);
            }
            TCLAP::ValuesConstraint<int> transfConst(transfVal);
            m_ttypeArg = new TCLAP::ValueArg<int>("", "ttype", transfStr, false, 0, &transfConst);
        
            // Create Quantization type description string
            std::vector<int> quantVal;
            std::string quantStr = "Quantization type : ";
            for(QuantizationType i = UNIFORM_QUANTIZATION; i <= NO_QUANTIZATION; i++)
            {
                quantStr += getQuantizationStringFromEnum(i).toStdString() + " : " + valToString(i);
                
                if(i == DEFAULT_QUANT_TYPE)
                {
                    quantStr += " (default)";
                }
                if(i < NO_QUANTIZATION)
                {
                    quantStr += ", ";
                }
                quantVal.push_back(i);
            }
            TCLAP::ValuesConstraint<int> quantConst(quantVal);
            m_qtypeArg = new TCLAP::ValueArg<int>("", "qtype", quantStr, false, 0, &quantConst);
        
            // Create Quantization number description string
            std::string qnumStr = "Number of digits for quantization : Uniform : rounding at 10^(qnum) position, ";
            qnumStr += "Non-uniform : keeps qnum most significant digits. (default ";
            qnumStr += valToString(DEFAULT_PRECISION_QUANT) + ")";
            m_qnumArg = new TCLAP::ValueArg<int>("", "qnum", qnumStr, false, 0, "int");
        
        
            // Create Encoding type description string
            std::vector<int> encVal;
            std::string encStr = "Encoding type : ";
            for(EncodingType i = HUFFMAN_ENCODING; i <= NO_ENCODING; i++)
            {
                encStr += getEncodingStringFromEnum(i).toStdString() + " : " + valToString(i);
                
                if(i == DEFAULT_ENCODING_TYPE)
                {
                    encStr += " (default)";
                }
                if(i < NO_ENCODING)
                {
                    encStr += ", ";
                }
                encVal.push_back(i);
            }
            TCLAP::ValuesConstraint<int> encConst(encVal);
            m_etypeArg = new TCLAP::ValueArg<int>("", "etype", encStr, false, 0, &encConst);
        
            // Add arguments to commandline
            cmd.add(m_ttypeArg);
            cmd.add(m_qtypeArg);
            cmd.add(m_qnumArg);
            cmd.add(m_etypeArg);
        #endif
        
        #if _SAVE_STATS
            m_statsArg = new TCLAP::ValueArg<std::string>("", "stats", "Statistics output text file", false, "", "string");
            cmd.add(m_statsArg);
        #endif

        cmd.parse(argc, argv);
        
    }
    catch( TCLAP::ArgException &e)
    {
        std::cerr << "error : " << e.error() << " for " << e.argId() << std::endl;
        return false;
    }
    return true;
}


/********************************************//**
\brief get parameters from commandline and validate them
\param parameters : structure to fill containing all parameters
\return boolean true if all arguments are valid, false otherwise
***********************************************/
bool Application::getAndValidateCommandLineArgs(Param& parameters)
{
    // Validate arguments
    QString tmp;
    parameters.inputPath = tmp.fromStdString(m_inputArg->getValue());
    std::ifstream input_file(m_inputArg->getValue().c_str());
    if(!input_file.good())
    {
        showMessage("Error : Invalid input file.");
        return false;
    }

    parameters.outputPath = tmp.fromStdString(m_outputArg->getValue());
    std::ifstream output_file(m_outputArg->getValue().c_str());
    if(output_file.good() && !m_forceArg->isSet())
    {
        showMessage("Error : Output file already exists. Use -f to overwrite output file.");
        return false;
    }
    
    parameters.compress = m_compArg->isSet();
    
    parameters.errorMax = DEFAULT_ERROR_MAX;
    if(m_errorArg->isSet())
    {
        if(m_errorArg->getValue() <= 0.0f)
        {
            showMessage("Error : Maximum error must greater than zero.");
            return false;
        }
        parameters.errorMax = m_errorArg->getValue();
    }
    
    parameters.encode = DEFAULT_ENCODE;
    parameters.transfType = DEFAULT_TRANSF_TYPE;
    parameters.quantizType = DEFAULT_QUANT_TYPE;
    parameters.precisionQuantiz = DEFAULT_PRECISION_QUANT;
    parameters.encodingType = DEFAULT_ENCODING_TYPE;
    parameters.statsOutputPath = "";
    
    #if !_USE_LIGHTWEIGHT
        if(m_ttypeArg->isSet())
        {
            parameters.transfType = (TransformationType)m_ttypeArg->getValue();
        }
        if(m_qtypeArg->isSet())
        {
            parameters.quantizType = (QuantizationType)m_qtypeArg->getValue();
        }
        
        if(m_qnumArg->isSet())
        {
            parameters.precisionQuantiz = m_qnumArg->getValue();
            if(parameters.quantizType == NO_QUANTIZATION && m_qnumArg->isSet())
            {
                showMessage("Error : A quantization precision is specified, but no quantization has been selected.");
                return false;
            }
            if(parameters.quantizType == NON_UNIFORM_QUANTIZATION && parameters.precisionQuantiz <= 0)
            {
                showMessage("Error : Number of coefficients for non-uniform quantization must be at least 1 or higher.");
                return false;
            }
        }
        
        if(m_etypeArg->isSet())
        {
            parameters.encodingType = (EncodingType)m_etypeArg->getValue();
        }
    #endif
    
    #if _SAVE_STATS
        if(m_outputArg->isSet())
        {
            parameters.statsOutputPath = tmp.fromStdString(m_statsArg->getValue());
        }
    #endif
    
    // Validate parameters
    QString errorTitle, errorMsg;
    if(!validateParameters(m_param, errorTitle, errorMsg))
    {
        showMessage(errorMsg);
        return false;
    }

    return true;
}


/********************************************//**
\brief start process from commandline
\param param : parameters to use in process
***********************************************/
void Application::startFromCommandLine(const Param& param)
{
    // Save parameters if needed
    #if _SAVE_STATS
        if(!param.statsOutputPath.isEmpty())
        {
            QFile myFile;
            myFile.setFileName(param.statsOutputPath);
            saveParameters(myFile, param);
        }
    #endif
    
    // Initialize fibers parameters
    m_pFibers->init(param.compress,
                    param.inputPath,
                    param.outputPath,
                    param.statsOutputPath,
                    param.errorMax,
                    param.transfType,
                    param.quantizType,
                    param.precisionQuantiz,
                    param.encodingType);
    
    // Start process
    m_timer.start();
    m_pThread->start();
    m_pThread->wait();
}


/********************************************//**
\brief close application and delete objects
***********************************************/
void Application::close()
{
    // Delete objets
    if(m_pThread != NULL)
    {
        m_pThread->exit(1);
        m_pThread->deleteLater();
    }
    
    // Exit from application
    QCoreApplication::exit(0);
}


/********************************************//**
\brief slot triggered when the the processing is finished
***********************************************/
void Application::processFinished()
{
    QString toPrint = "Time elapsed : " + QString::number(m_timer.elapsed() / 1000.0 / 60.0) + " min.\n";
    
    // Save time elapsed if neeeded
    #if _SAVE_STATS
        if(!m_param.statsOutputPath.isEmpty())
        {
            m_file.setFileName(m_param.statsOutputPath);
            QTextStream out(&m_file);
            
            if(m_file.open(QIODevice::WriteOnly | QIODevice::Text | QIODevice::Append))
            {
                out << toPrint;
            }
            m_file.close();
        }
    #endif
    
    // Print time elapsed
    showMessage(toPrint);

    // Close application
    close();
}
