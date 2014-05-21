
#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QButtonGroup>
#include <QComboBox>
#include <QDialogButtonBox>
#include <QElapsedTimer>
#include <QFile>
#include <QFileInfo>
#include <QGroupBox>
#include <QLabel>
#include <QLayout>
#include <QLineEdit>
#include <QMessageBox>
#include <QObject>
#include <QProgressBar>
#include <QRadioButton>
#include <QSpinBox>
#include <QString>
#include <QStringList>
#include <QtGui>
#include <QToolButton>
#include <QVariant>

#include "Fibers.h"
#include "Utils.h"


class MainWindow : public QMainWindow
{
    Q_OBJECT
    
public:
    MainWindow(QWidget *parent = 0);
    ~MainWindow();
    
public slots:
    void showMessage(const QString&);
    void processFinished();
    void stopProcess();
    
private slots:
    void updateUI();
    void resetButtonClicked();
    void startButtonClicked();
    void inputToolBtnClicked();
    void outputToolBtnClicked();
    void outputStatsToolBtnClicked();

signals:
    void stop();
    void cleanUp();
    void reset();
    
private:
    void connectSlots();
    void closeEvent(QCloseEvent * event);
    void createStatusBar();
    void createWindow();
    void disableControls(bool disable);
    void getLastDirectory();
    void startFromGUI();
    void updateLastDirectory();
    
    
    // Layouts
    QWidget* m_pWidget;
    QGridLayout* m_pLayout;
    QGridLayout* m_pOptionsLayout;
    QHBoxLayout* m_pAllOptionsLayout;
    QVBoxLayout* m_pFirstSectionLayout;
    
    // Input widgets and layouts
    QLabel* m_pInputText;
    QLineEdit* m_pInputTextBox;
    QPushButton* m_pInputToolBtn;
    QVBoxLayout* m_pInputLayout;
    QGridLayout* m_pInputWidgetsLayout;
    
    // Output widgets and layouts
    QLabel* m_pOutputText;
    QLineEdit* m_pOutputTextBox;
    QPushButton* m_pOutputToolBtn;
    QVBoxLayout* m_pOutputLayout;
    QGridLayout* m_pOutputWidgetsLayout;
    
    // Compression widgets and layouts
    QButtonGroup* m_pCompressBtnGroup;
    QGroupBox* m_pCompressionBox;
    QRadioButton* m_pCompressRadioBtn;
    QRadioButton* m_pDecompressRadioBtn;
    QVBoxLayout* m_pCompressLayout;
    QVBoxLayout* m_pCompressBoxLayout;
    
    // Maximum error widgets and layouts
    QGroupBox* m_pErrorBox;
    QLabel* m_pErrorMM;
    QLineEdit* m_pErrorTextBox;
    QHBoxLayout* m_pErrorLayout;
    
    #if !_USE_LIGHTWEIGHT
        // Transformation widgets and layouts
        QComboBox* m_pTransformationComboBox;
        QGroupBox* m_pTransformationBox;
        QVBoxLayout* m_pTransformationBoxLayout;
        
        // Quantization widgets and layouts
        QComboBox* m_pQuantizationComboBox;
        QGroupBox* m_pQuantizationBox;
        QSpinBox* m_pQuantizationSpinBox;
        QVBoxLayout* m_pQuantizationBoxLayout;
        
        // Encoding widgets and layouts
        QComboBox* m_pEncodingComboBox;
        QGroupBox* m_pEncodingBox;
        QVBoxLayout* m_pEncodingBoxLayout;
    #endif
    
    #if _SAVE_STATS
        // Output stats widgets and layouts
        QLabel* m_pOutputStatsText;
        QLineEdit* m_pOutputStatsTextBox;
        QPushButton* m_pOutputStatsToolBtn;
        QVBoxLayout* m_pOutputStatsLayout;
        QGridLayout* m_pOutputStatsWidgetsLayout;
    #endif
    
    // Dialog buttons widgets
    QDialogButtonBox* m_pDialogButtonBox;
    QPushButton* m_pResetBtn;
    QPushButton* m_pStartBtn;
    QHBoxLayout* m_pDialogLayout;
    
    // ProgressBar widget
    QProgressBar* m_pProgressBar;
    
    // StatusBar
    QLabel* m_pStatusLabel;
    QStatusBar* m_pStatusBar;
    QHBoxLayout* m_pStatusLayout;
    
    QString m_lastDirectory;
    Fibers* m_pFibers;
    QThread* m_pThread;
    Param m_param;
    
    QElapsedTimer m_timer;
    QFile m_file;
};

#endif // MAINWINDOW_H
