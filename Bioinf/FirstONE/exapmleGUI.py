import sys
from PyQt6.QtWidgets import QApplication, QWidget, QLabel, QLineEdit, QPushButton, QVBoxLayout, QMessageBox

class SimpleGUI(QWidget):
    def __init__(self):
        super().__init__()

        self.initUI()

    def initUI(self):
        # Tworzenie widgetów
        self.label1 = QLabel('Zmienna 1:', self)
        self.entry1 = QLineEdit(self)

        self.label2 = QLabel('Zmienna 2:', self)
        self.entry2 = QLineEdit(self)

        self.button = QPushButton('Wykonaj operację', self)
        self.button.clicked.connect(self.execute_operation)

        # Ustawienie layoutu
        vbox = QVBoxLayout()
        vbox.addWidget(self.label1)
        vbox.addWidget(self.entry1)
        vbox.addWidget(self.label2)
        vbox.addWidget(self.entry2)
        vbox.addWidget(self.button)

        self.setLayout(vbox)

        # Ustawienia okna
        self.setWindowTitle('Interfejs użytkownika')
        self.show()

    def execute_operation(self):
        variable1 = self.entry1.text()
        variable2 = self.entry2.text()
        result = f"Wprowadzone wartości to: {variable1} i {variable2}"
        QMessageBox.information(self, 'Wynik', result)

if __name__ == '__main__':
    app = QApplication(sys.argv)
    ex = SimpleGUI()
    sys.exit(app.exec_())
