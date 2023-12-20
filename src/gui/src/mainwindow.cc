#include "./ui_mainwindow.h"
#include "diffschemes.hpp"
#include "logger.hpp"
#include "mainwindow.h"
#include "nmlib.hpp"
#include <cmath>
#include <cstddef>
#include <exception>
#include <fstream>
#include <qcolor.h>
#include <qcustomplot.h>
#include <qvalidator.h>
#include <stdexcept>
#include <type_traits>

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
    QVector<double> x, v, u;

    for (int i = 0; i < N; i++) {
        x.push_back(res1.at(i).x);
        v.push_back(res1.at(i).v);
        u.push_back(res1.at(i).u);
    }

    this->ui->plot->addGraph();

    QPen pen;
    pen.setWidth(2);
    pen.setColor(col);
    this->ui->plot->graph(count_plot)->setPen(pen);
    this->ui->plot->addGraph();

    this->ui->plot->graph(count_plot)->addData(x, v);
    this->ui->plot->replot();
    count_plot++;

    if (func == 0) {
        pen.setWidth(2);
        pen.setColor(QColor(0x7f00ff));
        this->ui->plot->graph(count_plot)->setPen(pen);
        this->ui->plot->addGraph();

        this->ui->plot->graph(count_plot)->addData(x, u);
        this->ui->plot->replot();
        count_plot++;
    }
}

void MainWindow::on_button_clear_clicked() {
    for (int i = 0; i < count_plot; i++) {
        this->ui->plot->graph(i)->data()->clear();
    }
    count_plot = 0;
    this->ui->plot->replot();
    this->ui->plot->update();
}

void MainWindow::on_button_save_points_clicked() {}

void MainWindow::on_button_plot_from_file_clicked() {}

void MainWindow::on_exit_button_clicked() {
    QMessageBox::StandardButton exit = QMessageBox::question(this, " ", "you wanna exit?", QMessageBox::Yes | QMessageBox::No);
    if (exit == QMessageBox::Yes) {
        QApplication::quit();
    }
}

void MainWindow::on_getdata_buttom_clicked() {
   
    N = this->ui->n_te->text().toInt();

    net = Uniform1DNet{ 0.l, 1.l, N };

    res1 = Table{ N };

    switch (func) {
    case 0: {
        auto ca = [](const double& x, const double& y, const double& step) -> double {
            if (x == y)
                return 0.l;
            if (y < x)
                throw std::runtime_error("error while evaluating _ca");

            if (y <= 0.5) {
                return 9.l / 4.l * step * 1.l / (y - x);
            } else if (x < 0.5 && y > 0.5) {
                return step * 1.l / (4.l / 9.l * (0.5 - x) + 4.l * (y - 0.5));
            } else if (x >= 0.5) {
                return 1.l / 4.l * step * 1 / (y - x);
            }
        };

        auto cd = [](const double& x, const double& y, const double& step) -> double {
            if (x == y)
                return 0.l;
            if (y < x)
                throw std::runtime_error("error while evaluating _ca");

            if (x <= 1.l && y <= 1.l) {
                return (y - x) / step;
            }
        };

        auto cphi = [](const double& x, const double& y, const double& step) -> double {
            if (x == y)
                return 0.l;
            if (y < x)
                throw std::runtime_error("error while evaluating _ca");

            if (y <= 0.5) {
                return 0.l;
            } else if (x < 0.5 && y > 0.5) {
                return (y - 0.5l) / step;
            } else if (x >= 0.5) {
                return (y - x) / step;
            }
        };

        auto sol = [](const double& x) -> double {
            constexpr double c1 = 0.063933520772499621525817123839948954;
            constexpr double c2 = -c1;
            constexpr double c3 = -0.10143789939117362531986424993355064;
            constexpr double c4 = -1.8507344489243619669703272176533937;

            if (x <= 0.5) {
                return c1 * std::exp((2.l / 3) * x) + c2 * std::exp(-((2.l / 3) * x));
            } else if (x > 0.5) {
                return c3 * std::exp(2 * x) + c4 * std::exp(-2 * x) + 1.l;
            }
        };

        NMbalance method{ &net, &res1, std::move(ca), std::move(cd), std::move(cphi), 0.l, 0.l, std::move(sol) };
        method.eval();
        break;
    }
    case 1: {
        auto ca = [](const double& x, const double& y, const double& step) -> double {
            if (x == y)
                return 0.l;
            if (y < x)
                throw std::runtime_error("error while evaluating _ca");

            if (y <= 0.5) {
                return step * 1 / (1.l / (x + 1.l) - 1.l / (y + 1.l));
            } else if (y > 0.5 && x < 0.5) {
                return step * 1 / (1.l / (x + 1.l) - 1.l / y + 4.l / 3.l);
            } else if (x >= 0.5) {
                return step * 1 / (1.l / x - 1.l / y);
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
                constexpr double c = 2.l;
                return 1.l / step * (a * std::exp(-x) + b * std::exp(y) - c);
            } else if (x >= 0.5) {
                constexpr double a = 1 / std::exp(0.5l);
                return 1.l / step * a * (std::exp(y) - std::exp(x));
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
                return 1.l / step * (y - x);
            }
        };

        auto sol = [](const double& x) -> double {
            return 0.l;
        };

        NMbalance method{ &net, &res1, std::move(ca), std::move(cd), std::move(cphi), 0.l, 0.l, std::move(sol) };
        method.eval();
        break;
    }
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
    ui->tableWidget->setColumnCount(4);

    switch (func)
    {
    case 0:
        ui->tableWidget->setHorizontalHeaderLabels(QStringList()  << "xi" << "ui" << "vi" << "vi - ui");
        break;
    case 1:
        ui->tableWidget->setHorizontalHeaderLabels(QStringList()  << "xi" << "vi" << "v2i" << "vi(xi) - v2i(xi)");
        break;

    default:
        break;
    }

    for (int j = 0; j < N; j++) {
            QTableWidgetItem* item = new QTableWidgetItem(QString::number(res1.at(j).x));
            ui->tableWidget->setItem(j, 0, item);

            QTableWidgetItem* item1 = new QTableWidgetItem(QString::number(res1.at(j).v));
            ui->tableWidget->setItem(j, 1, item1);

            QTableWidgetItem* item2 = new QTableWidgetItem(QString::number(res1.at(j).u));
            ui->tableWidget->setItem(j, 2, item2);

            QTableWidgetItem* item3 = new QTableWidgetItem(QString::number(res1.at(j).eps));
            ui->tableWidget->setItem(j, 3, item3);

    }

}

void MainWindow::on_comboBox_activated(int index) {
    func = index;
}

void MainWindow::on_HelpButton_clicked() {
    form.show();
}
