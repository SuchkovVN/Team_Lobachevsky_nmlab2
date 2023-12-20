#include "./ui_mainwindow.h"
#include "diffschemes.hpp"
#include "logger.hpp"
#include "mainwindow.h"
#include "nmlib.hpp"
#include <cmath>
#include <cstddef>
#include <exception>
#include <fstream>
#include <qcustomplot.h>
#include <qvalidator.h>
#include <stdexcept> 

MainWindow::MainWindow(QWidget* parent) : QMainWindow(parent), ui(new Ui::MainWindow) {
    ui->setupUi(this);

    this->ui->plot->xAxis->setLabel("X");
    this->ui->plot->yAxis->setLabel("Y");

    this->ui->plot->xAxis->setRange(-10, 10);
    this->ui->plot->yAxis->setRange(-7, 7);

    this->ui->plot->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);

    func = 0;
}

MainWindow::~MainWindow() {
    delete ui;
}

void MainWindow::on_button_plot_clicked() {
    QVector<double> x, v;

    x_begin = this->ui->lineEdit_x_start->text().toDouble();
    x_end = this->ui->lineEdit_x_end->text().toDouble();
    h = this->ui->lineEdit_step->text().toDouble();
    precision = this->ui->lineEdit_precision->text().toDouble();

    for (int i = 0; i < res1.size(); i++) {
        x.push_back(res1.at(i).x);
        v.push_back(res1.at(i).v);
    }

    this->ui->plot->addGraph();

    QPen pen;
    pen.setWidth(2);
    pen.setColor(col);
    this->ui->plot->graph(count_plot)->setPen(pen);
    this->ui->plot->addGraph();
    pen.setWidth(2);
    pen.setColor(col);

    this->ui->plot->graph(count_plot)->addData(x, v);
    this->ui->plot->replot();
    count_plot++;
}

void MainWindow::on_button_clear_clicked() {
    for (int i = 0; i < count_plot; i++) {
        this->ui->plot->graph(i)->data()->clear();
    }
    count_plot = 0;
    this->ui->plot->replot();
    this->ui->plot->update();
}

void MainWindow::on_button_save_points_clicked() {
}

void MainWindow::on_button_plot_from_file_clicked() {}

void MainWindow::on_exit_button_clicked() {
    QMessageBox::StandardButton exit = QMessageBox::question(this, " ", "you wanna exit?", QMessageBox::Yes | QMessageBox::No);
    if (exit == QMessageBox::Yes) {
        QApplication::quit();
    }
}

void MainWindow::on_getdata_buttom_clicked() {
    x_begin = this->ui->lineEdit_x_start->text().toDouble();
    x_end = this->ui->lineEdit_x_end->text().toDouble();
    h = this->ui->lineEdit_step->text().toDouble();
    precision = this->ui->lineEdit_precision->text().toDouble();
    x_start = this->ui->lineEdit_start_x->text().toDouble();
    y_start = this->ui->lineEdit_start_y->text().toDouble();
    du = this->ui->lineEdit_last->text().toDouble();
    N = this->ui->lineEdit_n->text().toInt();
    A = this->ui->lineEdit_a->text().toDouble();
    B = this->ui->lineEdit_b->text().toDouble();
    C = this->ui->lineEdit_c->text().toDouble();

    net = Uniform1DNet{0.l, 1.l, N};

    res1 = Table{N};

    auto ca = [](const double& x, const double& y, const double& step) -> double {
        if (x == y)
            return 0.l;
        if (y < x)
            throw std::runtime_error("error while evaluating _ca");

        if (y <= 0.5) {
            return step * 1 / (1.l / (x + 1.l) - 1.l / (y + 1.l));
        } else if (y > 0.5 && x < 0.5) {
            return step * (1.l / (x + 1.l) - 1.l / y + 4.l / 3.l);
        } else if (x >= 0.5) {
            return step * (1.l / x - 1.l / y);
        }
    };

    auto cd = [](const double& x, const double& y, const double& step) -> double {
        if (x == y)
            return 0.l;
        if (y < x)
            throw std::runtime_error("error while evaluating _ca");

        if (y <= 0.5) {
            constexpr double a = std::exp(0.5l);
            return a / step * (std::exp(-x) - std::exp(-y));
        } else if (y > 0.5 && x < 0.5) {
            constexpr double a = std::exp(0.5l);
            constexpr double b = 1.l / a;
            constexpr double c = 2.l * std::exp(-0.5l) * std::exp(0.5l);
            return 1.l / step * (a * std::exp(-x) + b * std::exp(y) - c);
        } else if (x >= 0.5) {
            constexpr double a = std::exp(-0.5l);
            return step *  a * (std::exp(y) - std::exp(x));
        }
    };

    auto cphi = [](const double& x, const double& y, const double& step) -> double {
        if (x == y)
            return 0.l;
        if (y < x)
            throw std::runtime_error("error while evaluating _ca");
        
        constexpr double pi = 3.141592653589793238;
        constexpr double tpi = 1.l / pi;
        if (y <= 0.5) {
            return 1.l / step * tpi * (std::sin(pi * y) - std::sin(pi * x));
        } else if (y > 0.5 && x < 0.5) {
            constexpr double a = tpi - 0.5l;
            return 1.l / step * (y - std::sin(pi * x) * tpi + a);
        } else if (x >= 0.5) {
            return step * (y - x);
        } 
    };

    NMbalance method{&net, &res1, std::move(ca), std::move(cd), std::move(cphi), 0.l, 0.l};

    switch (func) {
    case 0:
        method.eval();
        break;
    case 1:
        
        break;
    default:
        break;
    }
}

//void MainWindow::on_Help_buttom_clicked(){};

void MainWindow::on_radioButton_blue_clicked(bool checked) {
    if (checked) {
        col = QColor(0, 0, 255);
    }
}

void MainWindow::on_radioButton_red_clicked(bool checked) {
    if (checked) {
        col = QColor(255, 0, 0);
    }
}

void MainWindow::on_radioButton_green_clicked(bool checked) {
    if (checked) {
        col = QColor(0, 255, 0);
    }
}

void MainWindow::on_radioButton_violet_clicked(bool checked) {
    if (checked) {
        col = QColor(105, 0, 198);
    }
}

void MainWindow::on_radioButton_mistake_clicked(bool checked) {
    LEC = checked;
}

void MainWindow::on_button_table_clicked() {
    ui->tableWidget->clear();

    ui->tableWidget->setRowCount(res1.size());
    ui->tableWidget->setColumnCount(14);

    std::cout << res1.size() << std::endl;
}

void MainWindow::on_comboBox_activated(int index) {
    func = index;
}

void MainWindow::on_HelpButton_clicked() {
    form.show();
}
