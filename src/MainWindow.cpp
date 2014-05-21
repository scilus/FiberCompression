
#include <QElapsedTimer>
#include "MainWindow.h"


/********************************************//**
\brief default constructor
\param parent : parent of the mainwindow
***********************************************/
MainWindow::MainWindow(QWidget *parent) :
QMainWindow(parent),
m_pFibers(NULL),
m_pThread(NULL)
{
    // Set central widget
    m_pWidget = new QWidget();
    setCentralWidget(m_pWidget);
    
    // Set title
    setWindowTitle(QString("FiberCompression"));
    
    // Set icon
    setWindowIcon(QIcon(":/icon.png"));
    setIconSize(QSize(32, 32));
    
    // Create window and controls
    createWindow();
    
    // Connect signals and slots
    connectSlots();
    
    // Get last directory
    getLastDirectory();
    
    // Update UI before starting
    updateUI();
}


/********************************************//**
\brief default destructor
***********************************************/
MainWindow::~MainWindow()
{
}


/********************************************//**
\brief connect signals and slots
***********************************************/
void MainWindow::connectSlots()
{
    connect(m_pInputToolBtn, SIGNAL(clicked()), this, SLOT(inputToolBtnClicked()));
    connect(m_pOutputToolBtn, SIGNAL(clicked()), this, SLOT(outputToolBtnClicked()));
    connect(m_pCompressRadioBtn, SIGNAL(clicked()), this, SLOT(updateUI()));
    connect(m_pDecompressRadioBtn, SIGNAL(clicked()), this, SLOT(updateUI()));
    #if !_USE_LIGHTWEIGHT
        connect(m_pQuantizationComboBox, SIGNAL(currentIndexChanged(int)), this, SLOT(updateUI()));
    #endif
    #if _SAVE_STATS
        connect(m_pOutputStatsToolBtn, SIGNAL(clicked()), this, SLOT(outputStatsToolBtnClicked()));
    #endif
    connect(m_pResetBtn, SIGNAL(clicked()), this, SLOT(resetButtonClicked()));
    connect(m_pStartBtn, SIGNAL(clicked()), this, SLOT(startButtonClicked()));
}


/********************************************//**
\brief create the window by adding controls
***********************************************/
void MainWindow::createWindow()
{
    // Input file widgets
    m_pInputText = new QLabel(QString("Input file (" + getStringValidExtensions() + ")"));
    m_pInputTextBox = new QLineEdit();
    m_pInputToolBtn = new QPushButton();
    m_pInputToolBtn->setText(QString("..."));
    
    // Input file layouts
    m_pInputLayout = new QVBoxLayout();
    m_pInputWidgetsLayout = new QGridLayout();
    m_pInputWidgetsLayout->addWidget(m_pInputTextBox, 1, 1);
    m_pInputWidgetsLayout->addWidget(m_pInputToolBtn, 1, 9);
    m_pInputLayout->addWidget(m_pInputText);
    m_pInputLayout->addLayout(m_pInputWidgetsLayout);
    
    // Output file widgets
    m_pOutputText = new QLabel(QString("Output file (" + getStringValidExtensions() + ")"));
    m_pOutputTextBox = new QLineEdit();
    m_pOutputToolBtn = new QPushButton();
    m_pOutputToolBtn->setText(QString("..."));
    
    // Output file layouts
    m_pOutputLayout = new QVBoxLayout();
    m_pOutputWidgetsLayout = new QGridLayout();
    m_pOutputWidgetsLayout->addWidget(m_pOutputTextBox, 1, 1);
    m_pOutputWidgetsLayout->addWidget(m_pOutputToolBtn, 1, 9);
    m_pOutputLayout->addWidget(m_pOutputText);
    m_pOutputLayout->addLayout(m_pOutputWidgetsLayout);
    
    // Compression / Decompression widgets
    m_pCompressionBox = new QGroupBox(QString("Mode"));
    m_pCompressBtnGroup = new QButtonGroup();
    m_pCompressRadioBtn = new QRadioButton(QString("Compress"));
    m_pCompressRadioBtn->setChecked(true);
    m_pDecompressRadioBtn = new QRadioButton(QString("Decompress"));
    m_pCompressBtnGroup->addButton(m_pCompressRadioBtn);
    m_pCompressBtnGroup->addButton(m_pDecompressRadioBtn);
    m_pCompressBoxLayout = new QVBoxLayout();
    m_pCompressBoxLayout->addWidget(m_pCompressRadioBtn);
    m_pCompressBoxLayout->addWidget(m_pDecompressRadioBtn);
    m_pCompressBoxLayout->setAlignment(m_pCompressBoxLayout, Qt::AlignVCenter);
    m_pCompressionBox->setLayout(m_pCompressBoxLayout);
    
    // Compression / Decompression layouts
    m_pCompressLayout = new QVBoxLayout();
    m_pCompressLayout->addWidget(m_pCompressionBox);
    
    // Maximum error widgets
    m_pErrorBox = new QGroupBox(QString("Maximum error"));
    m_pErrorTextBox = new QLineEdit();
    m_pErrorTextBox->setValidator(new QDoubleValidator());
    m_pErrorTextBox->setText(QString::number(DEFAULT_MAX_ERROR));
    m_pErrorMM = new QLabel(QString("mm"));
    
    // Maximum error layouts
    m_pErrorLayout = new QHBoxLayout();
    m_pErrorLayout->addWidget(m_pErrorTextBox);
    m_pErrorLayout->addWidget(m_pErrorMM);
    m_pErrorBox->setLayout(m_pErrorLayout);
    
    #if !_USE_LIGHTWEIGHT
        // Transformation widgets
        m_pTransformationBox = new QGroupBox(QString("Transformation"));
        m_pTransformationComboBox = new QComboBox();
        m_pTransformationComboBox->addItem(getTransformationStringFromEnum(NO_TRANSFORMATION));
        m_pTransformationComboBox->addItem(getTransformationStringFromEnum(DCT_TRANSFORMATION));
        m_pTransformationComboBox->addItem(getTransformationStringFromEnum(DB4_TRANSFORMATION));
        m_pTransformationComboBox->addItem(getTransformationStringFromEnum(DB6_TRANSFORMATION));
        m_pTransformationComboBox->addItem(getTransformationStringFromEnum(DB8_TRANSFORMATION));
        m_pTransformationComboBox->addItem(getTransformationStringFromEnum(BIOR_5_3_TRANSFORMATION));
        m_pTransformationComboBox->addItem(getTransformationStringFromEnum(BIOR_9_7_TRANSFORMATION));
        
        // Transformation layouts
        m_pTransformationBoxLayout = new QVBoxLayout();
        m_pTransformationBoxLayout->addWidget(m_pTransformationComboBox);
        m_pTransformationBox->setLayout(m_pTransformationBoxLayout);
        
        // Quantization widgets
        m_pQuantizationBox = new QGroupBox(QString("Quantization"));
        m_pQuantizationComboBox = new QComboBox();
        m_pQuantizationComboBox->addItem(QString("Uniform"));
        m_pQuantizationComboBox->addItem(QString("Non-Uniform"));
        m_pQuantizationComboBox->addItem(QString("No quantization"));
        m_pQuantizationSpinBox = new QSpinBox();
        
        // Quantization layouts
        m_pQuantizationBoxLayout = new QVBoxLayout();
        m_pQuantizationBoxLayout->addWidget(m_pQuantizationComboBox);
        m_pQuantizationBoxLayout->addWidget(m_pQuantizationSpinBox);
        m_pQuantizationBox->setLayout(m_pQuantizationBoxLayout);
        
        // Encoding widgets
        m_pEncodingBox = new QGroupBox(QString("Encoding"));
        m_pEncodingComboBox = new QComboBox();
        m_pEncodingComboBox->addItem(QString("Huffman"));
        m_pEncodingComboBox->addItem(QString("Arithmetic"));
        m_pEncodingComboBox->addItem(QString("No encoding"));
        
        // Encoding layouts
        m_pEncodingBoxLayout = new QVBoxLayout();
        m_pEncodingBoxLayout->addWidget(m_pEncodingComboBox);
        m_pEncodingBoxLayout->setAlignment(m_pEncodingComboBox, Qt::AlignTop);
        m_pEncodingBox->setLayout(m_pEncodingBoxLayout);
    #endif
    
    #if _SAVE_STATS
        // Input file widgets
        m_pOutputStatsText = new QLabel(QString("Output statistics file : "));
        m_pOutputStatsTextBox = new QLineEdit();
        m_pOutputStatsToolBtn = new QPushButton();
        m_pOutputStatsToolBtn->setText(QString("..."));
        
        // Input file layouts
        m_pOutputStatsLayout = new QVBoxLayout();
        m_pOutputStatsWidgetsLayout = new QGridLayout();
        m_pOutputStatsWidgetsLayout->addWidget(m_pOutputStatsTextBox, 1, 1);
        m_pOutputStatsWidgetsLayout->addWidget(m_pOutputStatsToolBtn, 1, 9);
        m_pOutputStatsLayout->addWidget(m_pOutputStatsText);
        m_pOutputStatsLayout->addLayout(m_pOutputStatsWidgetsLayout);
    #endif
    
    // Dialog buttons widgets
    m_pDialogButtonBox = new QDialogButtonBox();
    m_pResetBtn = new QPushButton(QString("Reset"));
    m_pStartBtn = new QPushButton(QString("Start"));
    m_pDialogButtonBox->addButton(m_pResetBtn, QDialogButtonBox::ResetRole);
    m_pDialogButtonBox->addButton(m_pStartBtn, QDialogButtonBox::AcceptRole);
    m_pDialogButtonBox->centerButtons();
    
    // Dialog buttons layout
    m_pDialogLayout = new QHBoxLayout();
    m_pDialogLayout->addWidget(m_pDialogButtonBox);
    m_pDialogLayout->setAlignment(m_pDialogButtonBox, Qt::AlignRight);
    
    // Progress bar widget
    m_pProgressBar = new QProgressBar();
    m_pProgressBar->hide();
    m_pProgressBar->setRange(0, 0);
    
    // Status Bar
    m_pStatusBar = new QStatusBar();
    m_pStatusBar->setSizeGripEnabled(false);
    
    m_pStatusLayout = new QHBoxLayout();
    m_pStatusLayout->addWidget(m_pStatusBar);
    m_pStatusLayout->addWidget(m_pProgressBar);
    
    // Options layout
    m_pOptionsLayout = new QGridLayout();
    m_pOptionsLayout->addWidget(m_pErrorBox, 1, 1);
    #if !_USE_LIGHTWEIGHT
        m_pOptionsLayout->addWidget(m_pTransformationBox, 1, 2);
        m_pOptionsLayout->addWidget(m_pQuantizationBox, 2, 1);
        m_pOptionsLayout->addWidget(m_pEncodingBox, 2, 2);
    #endif
    
    m_pAllOptionsLayout = new QHBoxLayout();
    m_pAllOptionsLayout->addLayout(m_pCompressLayout);
    m_pAllOptionsLayout->addLayout(m_pOptionsLayout);
    
    m_pFirstSectionLayout = new QVBoxLayout();
    m_pFirstSectionLayout->addLayout(m_pInputLayout);
    m_pFirstSectionLayout->addLayout(m_pOutputLayout);
    
    // Add layouts to main layout
    m_pLayout = new QGridLayout();
    m_pLayout->setContentsMargins(25, 20, 25, 20);
    m_pLayout->addLayout(m_pFirstSectionLayout, 1, 1);
    m_pLayout->addLayout(m_pAllOptionsLayout, 2, 1);
    
    #if _SAVE_STATS
        m_pLayout->addLayout(m_pOutputStatsLayout, 3, 1);
    #endif
    
    m_pLayout->addLayout(m_pDialogLayout, 4, 1);
    m_pLayout->addLayout(m_pStatusLayout, 5, 1);
    
    // Add widgets to central widget
    m_pWidget->setLayout(m_pLayout);
}


/********************************************//**
\brief update the UI when a control is modify
***********************************************/
void MainWindow::updateUI()
{
    // Disable all options if decompression is checked and enable all if compression is checked
    bool disableOptions = m_pDecompressRadioBtn->isChecked();
    disableControls(disableOptions);
    
    #if !_USE_LIGHTWEIGHT
        // Update min and max value of spin box according to quantization type
        if(m_pQuantizationComboBox->currentIndex() == 0)
        {
            m_pQuantizationSpinBox->setMinimum(MIN_UNIFORM);
            m_pQuantizationSpinBox->setMaximum(MAX_UNIFORM);
            m_pQuantizationSpinBox->setValue(DEFAULT_UNIFORM);
            m_pQuantizationSpinBox->setDisabled(false);
        }
        else if(m_pQuantizationComboBox->currentIndex() == 1)
        {
            m_pQuantizationSpinBox->setMinimum(MIN_NON_UNIFORM);
            m_pQuantizationSpinBox->setMaximum(MAX_NON_UNIFORM);
            m_pQuantizationSpinBox->setValue(DEFAULT_NON_UNIFORM);
            m_pQuantizationSpinBox->setDisabled(false);
        }
        else
        {
            m_pQuantizationSpinBox->setDisabled(true);
        }
    #endif
}


/********************************************//**
\brief update the text file containing the last directory used when loading files
***********************************************/
void MainWindow::updateLastDirectory()
{
    QFile file;
    file.setFileName("preferences.txt");
    
    if(file.open(QIODevice::WriteOnly | QIODevice::Text))
    {
        QTextStream out(&file);
        out << m_lastDirectory;
    }
    file.close();
}


/********************************************//**
\brief get last directory from the text file
***********************************************/
void MainWindow::getLastDirectory()
{
    QFile file;
    file.setFileName("preferences.txt");
    
    if(file.open(QIODevice::ReadOnly | QIODevice::Text))
    {
        QTextStream in(&file);
        m_lastDirectory = in.readLine();
    }
    file.close();
}


/********************************************//**
\brief event triggered when the input button is clicked
***********************************************/
void MainWindow::inputToolBtnClicked()
{
    QString inputPath = QFileDialog::getOpenFileName(this,
                                                     QString("Input File :"),
                                                     m_lastDirectory,
                                                     QString("Input file (*.trk *.tck *.zfib)"));
    m_pInputTextBox->setText(inputPath);
    m_lastDirectory = QFileInfo(inputPath).path();
    updateLastDirectory();
}


/********************************************//**
\brief event triggered when the output button is clicked
***********************************************/
void MainWindow::outputToolBtnClicked()
{
    QString outputPath = QFileDialog::getSaveFileName(this,
                                                      QString("Output File :"),
                                                      m_lastDirectory,
                                                      QString("Output file (*.trk *.tck *.zfib)"));
    m_pOutputTextBox->setText(outputPath);
    m_lastDirectory = QFileInfo(outputPath).path();
    updateLastDirectory();
}


/********************************************//**
\brief event triggered when the output stats button is clicked
***********************************************/
void MainWindow::outputStatsToolBtnClicked()
{
    #if _SAVE_STATS
        QString outputPath = QFileDialog::getSaveFileName(this,
                                                          QString("Output Statistics File :"),
                                                          m_lastDirectory,
                                                          QString("Output file (*.txt)"));
        m_pOutputStatsTextBox->setText(outputPath);
        m_lastDirectory = QFileInfo(outputPath).path();
        updateLastDirectory();
    #endif
}


/********************************************//**
\brief event triggered when the reset button is clicked
***********************************************/
void MainWindow::resetButtonClicked()
{
    // Reset input file
    m_pInputTextBox->setText(QString(""));
    
    // Reset output file
    m_pOutputTextBox->setText(QString(""));
    
    // Reset compression
    m_pCompressRadioBtn->setChecked(true);
    
    // Reset error max
    m_pErrorTextBox->setText(QString::number(DEFAULT_MAX_ERROR));
    
    #if !_USE_LIGHTWEIGHT
        // Reset transformation
        m_pTransformationComboBox->setCurrentIndex(0);
    
        // Reset quantization
        m_pQuantizationComboBox->setCurrentIndex(0);
    
        // Reset encoding
        m_pEncodingComboBox->setCurrentIndex(0);
    #endif
    
    #if _SAVE_STATS
        m_pOutputStatsTextBox->setText(QString(""));
    #endif
    
    disableControls(false);
}


/********************************************//**
\brief event triggered when the close button is clicked
***********************************************/
void MainWindow::closeEvent(QCloseEvent * event)
{
    event->ignore();
    if (QMessageBox::Yes == QMessageBox::question(this, "Close Confirmation?",
                                                  "Are you sure you want to exit?",
                                                  QMessageBox::Yes|QMessageBox::No))
    {
        emit stop();
        emit cleanUp();
        
        //QApplication::instance()->quit();
		event->accept();
    }
};


/********************************************//**
\brief disable (or enable) all controls
***********************************************/
void MainWindow::disableControls(bool disable)
{
    m_pInputTextBox->setDisabled(disable);
    m_pInputToolBtn->setDisabled(disable);
    m_pOutputTextBox->setDisabled(disable);
    m_pOutputToolBtn->setDisabled(disable);
    m_pErrorBox->setDisabled(disable);

    #if !_USE_LIGHTWEIGHT
        m_pTransformationBox->setDisabled(disable);
        m_pQuantizationBox->setDisabled(disable);
        m_pEncodingBox->setDisabled(disable);
    #endif
    
    #if _SAVE_STATS
        m_pOutputStatsTextBox->setDisabled(disable);
    #endif
}


/********************************************//**
\brief event triggered when the start (or stop) button is clicked
***********************************************/
void MainWindow::startButtonClicked()
{
    if(m_pStartBtn->text() == QString("Start"))
    {
        // Get parameters from GUI
        m_param.inputPath = m_pInputTextBox->text();
        m_param.outputPath = m_pOutputTextBox->text();
        m_param.compress = m_pCompressRadioBtn->isChecked();
        m_param.errorMax = m_pErrorTextBox->text().toFloat();
        
        m_param.transfType = DEFAULT_TRANSF_TYPE;
        m_param.quantizType = DEFAULT_QUANT_TYPE;
        m_param.precisionQuantiz = DEFAULT_PRECISION_QUANT;
        m_param.encodingType = DEFAULT_ENCODING_TYPE;
        m_param.encode = DEFAULT_ENCODE;
        m_param.statsOutputPath = "";
        
        #if !_USE_LIGHTWEIGHT
            m_param.transfType = getTransformationEnumFromString(m_pTransformationComboBox->currentText());
            m_param.quantizType = getQuantizationEnumFromString(m_pQuantizationComboBox->currentText());
            m_param.precisionQuantiz = m_pQuantizationSpinBox->value();
            m_param.encodingType = getEncodingEnumFromString(m_pEncodingComboBox->currentText());
            m_param.encode = (m_param.encodingType != NO_ENCODING);
        #endif
        
        #if _SAVE_STATS
            m_param.statsOutputPath = m_pOutputStatsTextBox->text();
        #endif
        
        // Validate parameters
        QString errorTitle, errorMsg;
        if(!validateParameters(m_param, errorTitle, errorMsg))
        {
            QMessageBox msgBox(QMessageBox::Critical,
                               errorTitle,
                               errorMsg,
                               QMessageBox::Ok,
                               NULL);
            msgBox.exec();
            return;

        }
        else
        {
            // Show progressbar
            m_pProgressBar->show();
            
            // Disable controls
            disableControls(true);
            
            // Change text to Stop
            m_pStartBtn->setText(QString("Stop"));
            
            // Start FiberCompression
            startFromGUI();
        }
    }
    else
    {
        emit stop();
    }
}


/********************************************//**
\brief slot triggered when the process has been stopped
***********************************************/
void MainWindow::stopProcess()
{
    emit reset();
    
    // Hide progress bar
    m_pProgressBar->hide();
    
    // Enable controls
    disableControls(false);
    
    // Show Message
    showMessage(QString("Stopped."));
    
    // Change text to Start
    m_pStartBtn->setText(QString("Start"));
}


/********************************************//**
\brief slot triggered when the status bar message must change
***********************************************/
void MainWindow::showMessage(const QString& message)
{
    QApplication::instance()->processEvents();
    m_pStatusBar->showMessage(message);
    QApplication::instance()->processEvents();
}


/********************************************//**
\brief start pipepine from the GUI (start button)
***********************************************/
void MainWindow::startFromGUI()
{
    showMessage(QString("Started."));
    m_pProgressBar->show();
    
    #if _SAVE_STATS
        if(!m_param.statsOutputPath.isEmpty())
        {
            m_file.setFileName(m_param.statsOutputPath);
            saveParameters(m_file, m_param);
        }
    #endif
    
    m_pFibers = NULL;
    m_pFibers = new Fibers();

    // Connect fibers slots
    connect(m_pFibers, SIGNAL(showMessage(const QString&)), this, SLOT(showMessage(const QString&)));
    connect(m_pFibers, SIGNAL(isStopped()), this, SLOT(stopProcess()));
    connect(m_pFibers, SIGNAL(finished()), this, SLOT(processFinished()));
    connect(m_pFibers, SIGNAL(finished()), m_pFibers, SLOT(deleteLater()));
    connect(this, SIGNAL(stop()), m_pFibers, SLOT(stop()));
    connect(this, SIGNAL(cleanUp()), m_pFibers, SLOT(cleanUp()));
    connect(this, SIGNAL(reset()), m_pFibers, SLOT(reset()));

    m_pFibers->init(m_param.compress,
                    m_param.inputPath,
                    m_param.outputPath,
                    m_param.statsOutputPath,
                    m_param.errorMax,
                    m_param.transfType,
                    m_param.quantizType,
                    m_param.precisionQuantiz,
                    m_param.encodingType);
    
    m_timer.start();
    
    // Create thread and start processing
    m_pThread = NULL;
    m_pThread = new QThread();
    m_pFibers->moveToThread(m_pThread);
    connect(m_pThread, SIGNAL(started()), m_pFibers, SLOT(process()));
    connect(m_pThread, SIGNAL(finished()), m_pThread, SLOT(deleteLater()));
    m_pThread->start();
}


/********************************************//**
\brief slot triggered when the process is finished
***********************************************/
void MainWindow::processFinished()
{
    // Hide progress bar
    m_pProgressBar->hide();
    
    // Stop button became start
    m_pStartBtn->setText(QString("Start"));
    
    // Enable controls
    disableControls(false);
    
    // Print time elapsed
    float time = m_timer.elapsed() / 1000.0 / 60.0;
    QString toPrint = "Time elapsed : " +
                      QString::number(time) +
                      " min.\n";
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
    
    showMessage(toPrint);
}